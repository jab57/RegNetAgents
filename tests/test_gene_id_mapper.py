#!/usr/bin/env python3
"""
Tests for GeneIDMapper symbol resolution.

Regression coverage for the synthetic-ID bug: gene_id_cache.pkl held real
Ensembl IDs in ensembl_to_symbol for every GREmLN network gene, but
symbol_to_ensembl held synthetic ENSG_CACHED_<symbol> placeholders for most of
them, so symbol queries reported real network genes as "not found".
"""

import os
import pickle
import sys

from regnetagents.gene_id_mapper import GeneIDMapper, SYNTHETIC_ID_PREFIX
from regnetagents_langgraph_mcp_server import get_workflow

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))
import build_network_cache  # noqa: E402


def test_no_placeholder_when_real_id_exists():
    """A symbol with a real Ensembl ID in the cache never resolves to a placeholder."""
    mapper = GeneIDMapper()
    real = {s.upper() for e, s in mapper.cache["ensembl_to_symbol"].items()
            if not e.startswith(SYNTHETIC_ID_PREFIX)}
    bad = [s for s in real
           if (mapper.symbol_to_ensembl(s) or "").startswith(SYNTHETIC_ID_PREFIX)]
    assert not bad, f"{len(bad)} symbols resolve to placeholders, e.g. {bad[:5]}"


def test_round_trip_for_real_ids():
    """symbol_to_ensembl(ensembl_to_symbol(id)) returns a real ID for the same symbol."""
    mapper = GeneIDMapper()
    for e, s in mapper.cache["ensembl_to_symbol"].items():
        if e.startswith(SYNTHETIC_ID_PREFIX):
            continue
        back = mapper.symbol_to_ensembl(s)
        assert back is not None and not back.startswith(SYNTHETIC_ID_PREFIX)
        assert (mapper.ensembl_to_symbol(back) or "").upper() == s.upper()


async def test_every_gremln_gene_resolves_by_symbol():
    """Every gene in the epithelial_cell network is reachable by its symbol."""
    workflow = await get_workflow()
    agent = workflow.modeling_agent
    mapper = agent.gene_mapper
    all_genes = set(agent.cache.network_indices["epithelial_cell"]["all_genes"])
    unresolved = []
    for e in all_genes:
        s = mapper.ensembl_to_symbol(e)
        if s is None:
            continue
        if mapper.symbol_to_ensembl(s) not in all_genes:
            unresolved.append(s)
    # The bug left ~9,400 unresolved. (ZNF724 has two in-network IDs; its symbol
    # resolves to one of them, which satisfies this check.)
    assert not unresolved, f"{len(unresolved)} unresolved, e.g. {unresolved[:10]}"


async def test_previously_missing_genes_now_found():
    """Genes previously reported as absent from the epithelial_cell network are found."""
    workflow = await get_workflow()
    agent = workflow.modeling_agent
    for gene in ["CDH1", "RNF43", "CCNE1", "ARID1A", "SNX7", "CDK5", "CCR7"]:
        result = agent.query_network("gene_neighbors", "epithelial_cell", gene=gene)
        assert not result.get("error"), f"{gene}: {result.get('message')}"
        assert result.get("gene") == gene


def test_cache_rebuild_repairs_placeholders(tmp_path):
    """build_network_cache.update_gene_id_cache replaces placeholders with known real IDs."""
    assert build_network_cache.SYNTHETIC_ID_PREFIX == SYNTHETIC_ID_PREFIX
    cache_path = tmp_path / "gene_id_cache.pkl"
    with open(cache_path, "wb") as f:
        pickle.dump({
            "symbol_to_ensembl": {"CDH1": f"{SYNTHETIC_ID_PREFIX}CDH1",
                                  "TP53": "ENSG00000141510"},
            "ensembl_to_symbol": {"ENSG00000039068": "CDH1",
                                  "ENSG00000141510": "TP53",
                                  f"{SYNTHETIC_ID_PREFIX}CDH1": "CDH1"},
        }, f)
    # Empty output dir: no network PKLs, so no MyGene.info lookups are attempted.
    build_network_cache.update_gene_id_cache(str(tmp_path / "no_networks"), str(cache_path))
    with open(cache_path, "rb") as f:
        s2e = pickle.load(f)["symbol_to_ensembl"]
    assert s2e["CDH1"] == "ENSG00000039068"
    assert s2e["TP53"] == "ENSG00000141510"
