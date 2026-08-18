"""
Tests for compare_network_contexts (issue #17).

Requires TCGA caches to be built:
    python scripts/build_tcga_cache.py --all
"""

import pytest
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from regnetagents_langgraph_workflow import RegNetAgentsWorkflow
from regnetagents.context_comparison import compare_network_contexts, _extract_symbols_with_weights


@pytest.fixture(scope="module")
def agent():
    workflow = RegNetAgentsWorkflow()
    return workflow.modeling_agent


def _tcga_available(agent) -> bool:
    return agent.tcga_cache is not None and bool(agent.tcga_cache.tcga_indices)


# ---------------------------------------------------------------------------
# Basic structure tests
# ---------------------------------------------------------------------------

def test_basic_comparison_myc_coad(agent):
    if not _tcga_available(agent):
        pytest.skip("TCGA cache not built — run scripts/build_tcga_cache.py --all")

    result = compare_network_contexts(agent, "MYC", "coad")

    assert result.get("error") is not True, result.get("message")
    assert result["gene"] == "MYC"
    assert result["population_averaged_context"] == "epithelial_cell"
    assert result["tumor_state_context"] == "tcga_coad"

    reg = result["regulators"]
    assert "conserved" in reg
    assert "tumor_state_only" in reg
    assert "population_averaged_only" in reg
    assert isinstance(reg["conserved"], list)
    assert isinstance(reg["conserved_count"], int)
    assert 0.0 <= reg["conserved_fraction"] <= 1.0

    tgt = result["targets"]
    assert "conserved_count" in tgt
    assert 0.0 <= tgt["conserved_fraction"] <= 1.0

    interp = result["interpretation"]
    assert interp["regulatory_rewiring"] in ("low", "moderate", "high")
    assert "conserved_fraction_regulators" in interp
    assert "tumor_specific_regulator_count" in interp


def test_conserved_fraction_valid_range(agent):
    if not _tcga_available(agent):
        pytest.skip("TCGA cache not built")

    result = compare_network_contexts(agent, "MYC", "coad")
    assert result.get("error") is not True
    frac = result["regulators"]["conserved_fraction"]
    assert 0.0 <= frac <= 1.0


def test_tp53_brca_returns_valid_structure(agent):
    """TP53 is present in both GREmLN epithelial and TCGA BRCA — comparison should succeed."""
    if not _tcga_available(agent):
        pytest.skip("TCGA cache not built")

    result = compare_network_contexts(agent, "TP53", "brca")
    assert result.get("error") is not True, result.get("message")
    assert result["gene"] == "TP53"
    assert "conserved_count" in result["regulators"]
    assert result["interpretation"]["regulatory_rewiring"] in ("low", "moderate", "high")


# ---------------------------------------------------------------------------
# Rewiring classification thresholds
# ---------------------------------------------------------------------------

def test_rewiring_classification_low():
    """Fraction >= 0.6 → low."""
    from regnetagents.context_comparison import compare_network_contexts as _cmp

    class _FakeAgent:
        def query_network(self, query_type, gene=None, cell_type=None,
                          network_source=None, tcga_network=None, **kw):
            # Both contexts return the same 10 regulators → fraction = 1.0
            return {
                "error": False,
                "regulators": [{"gene": f"TF{i}"} for i in range(10)],
                "targets":    [],
            }

    result = _cmp(_FakeAgent(), "GENE", "brca")
    assert result["interpretation"]["regulatory_rewiring"] == "low"


def test_rewiring_classification_high():
    """Fraction < 0.3 → high."""
    from regnetagents.context_comparison import compare_network_contexts as _cmp

    class _FakeAgent:
        def __init__(self):
            self._call = 0

        def query_network(self, query_type, gene=None, cell_type=None,
                          network_source=None, tcga_network=None, **kw):
            self._call += 1
            # First call (pop-averaged): TF0–TF9; second call (TCGA): TF10–TF19
            start = 0 if self._call == 1 else 10
            return {
                "error": False,
                "regulators": [{"gene": f"TF{i}"} for i in range(start, start + 10)],
                "targets":    [],
            }

    result = _cmp(_FakeAgent(), "GENE", "brca")
    assert result["interpretation"]["regulatory_rewiring"] == "high"


# ---------------------------------------------------------------------------
# Error handling
# ---------------------------------------------------------------------------

def test_graceful_error_invalid_cancer_type(agent):
    """Unknown cancer type → error dict, no exception."""
    result = compare_network_contexts(agent, "TP53", "gbm")
    assert result.get("error") is True
    assert "message" in result


def test_graceful_error_gene_not_in_network(agent):
    """Gene absent from GREmLN → error dict, no exception."""
    result = compare_network_contexts(agent, "NOTAREALGENE999", "brca")
    assert result.get("error") is True


# ---------------------------------------------------------------------------
# MI edge-weight field (tumor_state_only_weights)
# ---------------------------------------------------------------------------

def test_extract_symbols_with_weights_dict_input():
    entries = [{"gene": "TF1", "likelihood": 0.5}, {"gene": "TF2"}]
    symbols, weights = _extract_symbols_with_weights(entries)
    assert symbols == {"TF1", "TF2"}
    assert weights == {"TF1": 0.5, "TF2": None}


def test_extract_symbols_with_weights_legacy_input():
    """Legacy plain-string entries (no dict/weight info) → empty weights dict, no error."""
    symbols, weights = _extract_symbols_with_weights(["TF1", "TF2"])
    assert symbols == {"TF1", "TF2"}
    assert weights == {}


def test_extract_symbols_with_weights_empty_input():
    symbols, weights = _extract_symbols_with_weights([])
    assert symbols == set()
    assert weights == {}


def test_tumor_state_only_weights_keys_match_tumor_state_only():
    """tumor_state_only_weights must have exactly the same keys as tumor_state_only — no more, no fewer."""
    class _FakeAgent:
        def __init__(self):
            self._call = 0

        def query_network(self, query_type, gene=None, cell_type=None,
                          network_source=None, tcga_network=None, **kw):
            self._call += 1
            if self._call == 1:
                # population-averaged: TF0-4
                return {"error": False, "regulators": [{"gene": f"TF{i}"} for i in range(5)], "targets": []}
            # tumor-state: TF3-9 (TF5-9 are tumor-only), with MI weights
            return {
                "error": False,
                "regulators": [{"gene": f"TF{i}", "likelihood": round(i * 0.1, 2)} for i in range(3, 10)],
                "targets": [],
            }

    result = compare_network_contexts(_FakeAgent(), "GENE", "brca")
    reg = result["regulators"]
    assert set(reg["tumor_state_only_weights"].keys()) == set(reg["tumor_state_only"])
    assert reg["tumor_state_only_weights"]["TF9"] == 0.9


def test_tumor_state_only_weights_legacy_string_entries_no_error():
    """Legacy plain-string regulator entries must not raise and must yield None weights."""
    class _FakeAgent:
        def __init__(self):
            self._call = 0

        def query_network(self, query_type, gene=None, cell_type=None,
                          network_source=None, tcga_network=None, **kw):
            self._call += 1
            if self._call == 1:
                return {"error": False, "regulators": [], "targets": []}
            return {"error": False, "regulators": ["TF1", "TF2"], "targets": []}

    result = compare_network_contexts(_FakeAgent(), "GENE", "brca")
    reg = result["regulators"]
    assert reg["tumor_state_only"] == ["TF1", "TF2"]
    assert reg["tumor_state_only_weights"] == {"TF1": None, "TF2": None}


def test_tumor_state_only_weights_present_on_real_query(agent):
    """Integration check against the real TCGA/GREmLN caches."""
    if not _tcga_available(agent):
        pytest.skip("TCGA cache not built")

    result = compare_network_contexts(agent, "MYC", "coad")
    assert result.get("error") is not True
    reg = result["regulators"]
    assert "tumor_state_only_weights" in reg
    assert set(reg["tumor_state_only_weights"].keys()) == set(reg["tumor_state_only"])


# ---------------------------------------------------------------------------
# Regression: existing tools unaffected
# ---------------------------------------------------------------------------

def test_query_network_unaffected(agent):
    """query_network still works correctly after adding compare_network_contexts."""
    result = agent.query_network(
        query_type="gene_neighbors",
        gene="TP53",
        cell_type="epithelial_cell",
        network_source="cell_type",
    )
    assert result.get("error") is not True
