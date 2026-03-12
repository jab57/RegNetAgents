#!/usr/bin/env python3
"""
Tests for TCGA ARACNe network integration (issue #16).

Test strategy
-------------
Unit tests (always run):
    Use a minimal in-memory fake network injected directly into the agent's
    tcga_cache, so they pass without any downloaded TCGA data.

Integration tests (skipped when TCGA PKLs are absent):
    Test against the real BRCA cache if it is present on disk.  Run these
    after executing `python scripts/build_tcga_cache.py --cancer-type brca`.

All tests follow the same async + get_workflow() pattern used elsewhere in
this test suite.
"""

import os
import pytest
from unittest.mock import MagicMock
from regnetagents_langgraph_mcp_server import get_workflow

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_BRCA_PKL = os.path.join("models", "networks", "tcga", "brca", "network_index.pkl")
_BRCA_AVAILABLE = os.path.exists(_BRCA_PKL)

# Minimal synthetic network injected for unit tests
#   TP53 → (CDKN1A, MDM2, BAX)  with Likelihood / MoA
#   MDM2 → TP53  (feedback loop)
_FAKE_NETWORK = {
    "regulator_targets":     {"TP53": ["CDKN1A", "MDM2", "BAX"],
                               "MDM2": ["TP53"]},
    "target_regulators":     {"CDKN1A": ["TP53"], "MDM2": ["TP53"],
                               "BAX": ["TP53"], "TP53": ["MDM2"]},
    "all_genes":             ["BAX", "CDKN1A", "MDM2", "TP53"],
    "regulator_target_mi":   {"TP53": {"CDKN1A": 0.30, "MDM2": 0.25, "BAX": 0.20},
                               "MDM2": {"TP53": 0.15}},
    "regulator_target_moa":  {"TP53": {"CDKN1A": 1.0, "MDM2": -1.0, "BAX": 1.0},
                               "MDM2": {"TP53": -1.0}},
    "regulator_target_count": {"TP53": {"CDKN1A": 0, "MDM2": 0, "BAX": 0},
                                "MDM2": {"TP53": 0}},
    "pagerank_normalized":   {"TP53": 1.0, "MDM2": 0.6, "CDKN1A": 0.3, "BAX": 0.3},
    "num_genes":             4,
    "num_edges":             4,
    "num_regulons":          2,
    "id_type":               "symbol",
    "source":                "tcga_aracne",
    "cancer_type":           "fake",
    "cache_version":         4,
}


def _inject_fake_network(agent, cancer_type: str = "fake") -> None:
    """Inject _FAKE_NETWORK into agent.tcga_cache for unit testing."""
    if agent.tcga_cache is None:
        from regnetagents_langgraph_workflow import TCGANetworkCache
        agent.tcga_cache = TCGANetworkCache()
    agent.tcga_cache.tcga_indices[cancer_type] = _FAKE_NETWORK


# ---------------------------------------------------------------------------
# Unit tests — registry and loader
# ---------------------------------------------------------------------------

def test_tcga_registry_has_8_types():
    from regnetagents.tcga_registry import TCGA_NETWORK_REGISTRY, TCGA_CANCER_TYPES
    assert len(TCGA_NETWORK_REGISTRY) == 8
    assert sorted(TCGA_CANCER_TYPES) == ["brca", "coad", "hnsc", "luad",
                                          "lusc", "ov", "prad", "ucec"]


def test_tcga_registry_paths_are_correct_format():
    from regnetagents.tcga_registry import TCGA_NETWORK_REGISTRY
    for ct, entry in TCGA_NETWORK_REGISTRY.items():
        assert "label" in entry
        assert "csv" in entry
        assert entry["csv"].endswith(f"tcga/{ct}/network.csv")


def test_load_network_raises_on_unknown_type():
    from regnetagents.network_loader import load_network
    with pytest.raises(KeyError, match="Unknown TCGA cancer type"):
        load_network("fake_cancer")


def test_load_network_raises_filenotfound_when_csv_absent():
    from regnetagents.network_loader import load_network
    with pytest.raises(FileNotFoundError):
        load_network("brca")  # CSV not downloaded in CI


def test_tcga_cancer_type_enum():
    from regnetagents_langgraph_workflow import TCGACancerType
    values = {t.value for t in TCGACancerType}
    assert values == {"brca", "coad", "hnsc", "luad", "lusc", "ov", "prad", "ucec"}


# ---------------------------------------------------------------------------
# Unit tests — query_network (fake network injected)
# ---------------------------------------------------------------------------

async def test_tcga_query_network_stats():
    """network_stats returns correct counts from fake network."""
    workflow = await get_workflow()
    agent = workflow.modeling_agent
    _inject_fake_network(agent)

    result = agent.query_network("network_stats",
                                 network_source="tcga", tcga_network="fake")
    assert result.get("query_type") == "network_stats"
    assert result.get("network_source") == "tcga"
    assert result.get("tcga_network") == "fake"
    assert result["num_genes"] == 4
    assert result["num_edges"] == 4
    assert result["num_regulons"] == 2
    assert result["has_moa"] is True


async def test_tcga_query_network_top_regulators():
    """top_regulators returns TP53 as highest out-degree node."""
    workflow = await get_workflow()
    agent = workflow.modeling_agent
    _inject_fake_network(agent)

    result = agent.query_network("top_regulators",
                                 network_source="tcga", tcga_network="fake", top_n=2)
    assert result.get("query_type") == "top_regulators"
    assert isinstance(result["results"], list)
    assert len(result["results"]) == 2
    assert result["results"][0]["gene"] == "TP53"
    assert result["results"][0]["num_targets"] == 3


async def test_tcga_query_network_top_targets():
    """top_targets returns TP53 (regulated by MDM2) as highest in-degree."""
    workflow = await get_workflow()
    agent = workflow.modeling_agent
    _inject_fake_network(agent)

    result = agent.query_network("top_targets",
                                 network_source="tcga", tcga_network="fake")
    assert result.get("query_type") == "top_targets"
    assert any(r["gene"] == "MDM2" for r in result["results"])


async def test_tcga_query_network_gene_neighbors():
    """gene_neighbors for TP53 returns targets with likelihood and moa fields."""
    workflow = await get_workflow()
    agent = workflow.modeling_agent
    _inject_fake_network(agent)

    result = agent.query_network("gene_neighbors", gene="TP53",
                                 network_source="tcga", tcga_network="fake")
    assert result.get("query_type") == "gene_neighbors"
    assert result.get("gene") == "TP53"
    assert result.get("network_source") == "tcga"
    assert result["num_targets"] == 3

    for entry in result["targets"]:
        assert "gene" in entry
        assert "likelihood" in entry
        assert "moa" in entry

    # CDKN1A should be activating (+1)
    cdkn1a = next(t for t in result["targets"] if t["gene"] == "CDKN1A")
    assert cdkn1a["moa"] == 1.0
    # MDM2 should be repressive (-1)
    mdm2 = next(t for t in result["targets"] if t["gene"] == "MDM2")
    assert mdm2["moa"] == -1.0


async def test_tcga_gene_not_found():
    """gene_neighbors for unknown gene returns error dict."""
    workflow = await get_workflow()
    agent = workflow.modeling_agent
    _inject_fake_network(agent)

    result = agent.query_network("gene_neighbors", gene="FAKEGENE",
                                 network_source="tcga", tcga_network="fake")
    assert result.get("error") is True
    assert "not found" in result.get("message", "").lower()


async def test_tcga_query_missing_network():
    """Querying a cancer type not loaded returns informative error."""
    workflow = await get_workflow()
    agent = workflow.modeling_agent
    _inject_fake_network(agent)

    result = agent.query_network("network_stats",
                                 network_source="tcga", tcga_network="ov")
    assert result.get("error") is True
    assert "ov" in result.get("message", "")


async def test_tcga_query_missing_tcga_network_param():
    """Omitting tcga_network when network_source='tcga' returns error."""
    workflow = await get_workflow()
    agent = workflow.modeling_agent
    _inject_fake_network(agent)

    result = agent.query_network("network_stats", network_source="tcga")
    assert result.get("error") is True
    assert "tcga_network" in result.get("message", "")


async def test_tcga_confidence_filter_medium():
    """confidence_level='medium' filters out edges with likelihood <= 0.05."""
    workflow = await get_workflow()
    agent = workflow.modeling_agent
    _inject_fake_network(agent)

    # All fake edges have likelihood > 0.05 so medium filter keeps all three
    result = agent.query_network("gene_neighbors", gene="TP53",
                                 network_source="tcga", tcga_network="fake",
                                 confidence_level="medium")
    assert result["num_targets"] == 3


async def test_tcga_confidence_filter_high():
    """confidence_level='high' keeps only edges with likelihood > 0.1."""
    workflow = await get_workflow()
    agent = workflow.modeling_agent
    _inject_fake_network(agent)

    # TP53→CDKN1A (0.30), →MDM2 (0.25), →BAX (0.20) all > 0.1
    result = agent.query_network("gene_neighbors", gene="TP53",
                                 network_source="tcga", tcga_network="fake",
                                 confidence_level="high")
    assert result["num_targets"] == 3


# ---------------------------------------------------------------------------
# Unit tests — find_master_regulators (fake network injected)
# ---------------------------------------------------------------------------

async def test_tcga_find_master_regulators_basic():
    """find_master_regulators on TCGA returns TP53 for known target set."""
    workflow = await get_workflow()
    agent = workflow.modeling_agent
    _inject_fake_network(agent)

    result = agent.find_master_regulators(
        gene_set=["CDKN1A", "MDM2", "BAX"],
        network_source="tcga",
        tcga_network="fake",
        top_n=5,
    )
    assert "master_regulators" in result, f"Unexpected: {result}"
    assert len(result["master_regulators"]) > 0
    top = result["master_regulators"][0]
    assert top["gene"] == "TP53"


async def test_tcga_find_master_regulators_structure():
    """Each entry has the expected fields (no ensembl_id for symbol-keyed network)."""
    workflow = await get_workflow()
    agent = workflow.modeling_agent
    _inject_fake_network(agent)

    result = agent.find_master_regulators(
        gene_set=["CDKN1A", "MDM2", "BAX"],
        network_source="tcga",
        tcga_network="fake",
        top_n=3,
    )
    for entry in result["master_regulators"]:
        assert "rank" in entry
        assert "gene" in entry
        assert "regulon_size" in entry
        assert "overlap_count" in entry
        assert "enrichment_score" in entry
        assert "p_value" in entry
        assert "overlapping_genes" in entry
        assert entry["overlap_count"] > 0
        assert 0.0 <= entry["p_value"] <= 1.0


async def test_tcga_find_master_regulators_query_summary():
    """query_summary reflects tcga context correctly."""
    workflow = await get_workflow()
    agent = workflow.modeling_agent
    _inject_fake_network(agent)

    result = agent.find_master_regulators(
        gene_set=["CDKN1A", "MDM2", "BAX"],
        network_source="tcga",
        tcga_network="fake",
    )
    summary = result["query_summary"]
    assert summary["network_source"] == "tcga"
    assert summary["tcga_network"] == "fake"
    assert summary["gene_set_size"] == 3
    assert summary["genes_found_in_network"] == 3
    assert summary["network_size"] == 4


async def test_tcga_find_master_regulators_unknown_genes():
    """All-unknown gene set returns error response."""
    workflow = await get_workflow()
    agent = workflow.modeling_agent
    _inject_fake_network(agent)

    result = agent.find_master_regulators(
        gene_set=["FAKEGENE1", "FAKEGENE2"],
        network_source="tcga",
        tcga_network="fake",
    )
    assert result.get("error") is True


# ---------------------------------------------------------------------------
# GREmLN path unaffected by TCGA changes
# ---------------------------------------------------------------------------

async def test_gremln_query_network_still_works():
    """Default network_source='cell_type' route is unaffected by TCGA additions."""
    workflow = await get_workflow()
    agent = workflow.modeling_agent
    result = agent.query_network("network_stats", cell_type="epithelial_cell")
    assert result.get("query_type") == "network_stats"
    assert result["num_genes"] > 0
    assert "network_source" not in result  # GREmLN result has no network_source key


async def test_gremln_find_master_regulators_still_works():
    """Default network_source='cell_type' for find_master_regulators is unaffected."""
    workflow = await get_workflow()
    agent = workflow.modeling_agent
    result = agent.find_master_regulators(
        ["CDKN1A", "MDM2", "BAX"], cell_type="epithelial_cell", top_n=5
    )
    assert "master_regulators" in result
    summary = result["query_summary"]
    assert summary["cell_type"] == "epithelial_cell"
    assert "network_source" not in summary  # GREmLN path doesn't add network_source


# ---------------------------------------------------------------------------
# Integration tests — real BRCA cache (skipped if not present)
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not _BRCA_AVAILABLE, reason="TCGA BRCA cache not built")
async def test_tcga_brca_cache_loads():
    """Real BRCA cache loads and has expected structure."""
    workflow = await get_workflow()
    agent = workflow.modeling_agent
    data = agent.tcga_cache.tcga_indices.get("brca", {})
    assert data, "BRCA cache is empty"
    assert data.get("id_type") == "symbol"
    assert data.get("source") == "tcga_aracne"
    assert data.get("cancer_type") == "brca"
    assert len(data.get("all_genes", [])) > 1000


@pytest.mark.skipif(not _BRCA_AVAILABLE, reason="TCGA BRCA cache not built")
async def test_tcga_brca_network_stats():
    """BRCA network_stats returns a non-trivial network."""
    workflow = await get_workflow()
    result = workflow.modeling_agent.query_network(
        "network_stats", network_source="tcga", tcga_network="brca"
    )
    assert result.get("query_type") == "network_stats"
    assert result["num_genes"] > 1000
    assert result["num_edges"] > 10000
    assert result["has_moa"] is True


@pytest.mark.skipif(not _BRCA_AVAILABLE, reason="TCGA BRCA cache not built")
async def test_tcga_brca_tp53_neighbors():
    """TP53 should have targets and regulators in the BRCA network."""
    workflow = await get_workflow()
    result = workflow.modeling_agent.query_network(
        "gene_neighbors", gene="TP53",
        network_source="tcga", tcga_network="brca"
    )
    assert not result.get("error"), result.get("message")
    assert result["num_targets"] > 0 or result["num_regulators"] > 0
    # TCGA edges carry moa field
    for t in result["targets"][:5]:
        assert "moa" in t
        assert "likelihood" in t


@pytest.mark.skipif(not _BRCA_AVAILABLE, reason="TCGA BRCA cache not built")
async def test_tcga_brca_master_regulators():
    """Known TP53 targets produce a ranked master regulator list in BRCA."""
    tp53_targets = ["CDKN1A", "MDM2", "BAX", "GADD45A", "PUMA",
                    "BBC3", "CCNG1", "SESN2", "DDB2", "POLK"]
    workflow = await get_workflow()
    result = workflow.modeling_agent.find_master_regulators(
        gene_set=tp53_targets,
        network_source="tcga",
        tcga_network="brca",
        top_n=10,
    )
    assert "master_regulators" in result, f"Unexpected: {result}"
    assert len(result["master_regulators"]) > 0
    p_values = [r["p_value"] for r in result["master_regulators"]]
    assert p_values == sorted(p_values), "Results not sorted by p-value"
