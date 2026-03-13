"""
context_comparison.py — Compare regulatory wiring across network contexts.

Computes the overlap between population-averaged (GREmLN) and tumor-state (TCGA)
regulatory programs for a given gene. All logic is rule-based — no LLM required.
"""

from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from regnetagents_langgraph_workflow import RegNetAgentsModelingAgent


def compare_network_contexts(
    agent: "RegNetAgentsModelingAgent",
    gene: str,
    cancer_type: str,
    cell_type: str = "epithelial_cell",
) -> dict:
    """
    Compare regulatory wiring for a gene across population-averaged and tumor-state networks.

    Args:
        agent:       RegNetAgentsModelingAgent instance (provides query_network).
        gene:        Gene symbol (e.g. "MYC").
        cancer_type: TCGA cancer type key (e.g. "coad", "brca").
        cell_type:   GREmLN reference context (default: "epithelial_cell").

    Returns:
        Structured comparison dict with conserved/context-specific regulators and targets,
        overlap statistics, and a rule-based regulatory rewiring classification.
    """
    gene_upper = gene.strip().upper()

    # --- Query population-averaged (GREmLN) network --------------------------
    pop_result = agent.query_network(
        query_type="gene_neighbors",
        gene=gene_upper,
        cell_type=cell_type,
        network_source="cell_type",
    )

    if pop_result.get("error"):
        return {
            "error": True,
            "gene": gene_upper,
            "message": f"Population-averaged query failed: {pop_result.get('message', 'unknown error')}",
        }

    # --- Query tumor-state (TCGA) network ------------------------------------
    tumor_result = agent.query_network(
        query_type="gene_neighbors",
        gene=gene_upper,
        network_source="tcga",
        tcga_network=cancer_type,
    )

    if tumor_result.get("error"):
        return {
            "error": True,
            "gene": gene_upper,
            "message": f"Tumor-state query failed: {tumor_result.get('message', 'unknown error')}",
        }

    # --- Extract regulator and target sets -----------------------------------
    pop_neighbors   = pop_result.get("gene_neighbors", {})
    tumor_neighbors = tumor_result.get("gene_neighbors", {})

    # GREmLN gene_neighbors returns lists of dicts with "gene" key
    def _extract_symbols_from_gremln(entries) -> set:
        if not entries:
            return set()
        if isinstance(entries[0], dict):
            return {e["gene"] for e in entries if "gene" in e}
        return set(entries)

    # TCGA gene_neighbors returns lists of dicts with "gene" key (same shape)
    def _extract_symbols_from_tcga(entries) -> set:
        if not entries:
            return set()
        if isinstance(entries[0], dict):
            return {e["gene"] for e in entries if "gene" in e}
        return set(entries)

    pop_regulators   = _extract_symbols_from_gremln(pop_neighbors.get("regulators", []))
    pop_targets      = _extract_symbols_from_gremln(pop_neighbors.get("targets", []))
    tumor_regulators = _extract_symbols_from_tcga(tumor_neighbors.get("regulators", []))
    tumor_targets    = _extract_symbols_from_tcga(tumor_neighbors.get("targets", []))

    # --- Compute overlaps ----------------------------------------------------
    reg_conserved    = sorted(pop_regulators & tumor_regulators)
    reg_pop_only     = sorted(pop_regulators - tumor_regulators)
    reg_tumor_only   = sorted(tumor_regulators - pop_regulators)

    tgt_conserved    = sorted(pop_targets & tumor_targets)
    tgt_pop_only     = sorted(pop_targets - tumor_targets)
    tgt_tumor_only   = sorted(tumor_targets - pop_targets)

    # Conserved fraction relative to the union (Jaccard-style) for regulators
    reg_union = pop_regulators | tumor_regulators
    reg_conserved_fraction = (
        round(len(reg_conserved) / len(reg_union), 4) if reg_union else 0.0
    )

    tgt_union = pop_targets | tumor_targets
    tgt_conserved_fraction = (
        round(len(tgt_conserved) / len(tgt_union), 4) if tgt_union else 0.0
    )

    # --- Regulatory rewiring classification (rule-based) ---------------------
    if reg_conserved_fraction >= 0.6:
        rewiring = "low"
    elif reg_conserved_fraction >= 0.3:
        rewiring = "moderate"
    else:
        rewiring = "high"

    return {
        "gene": gene_upper,
        "population_averaged_context": cell_type,
        "tumor_state_context": f"tcga_{cancer_type}",
        "regulators": {
            "population_averaged_total": len(pop_regulators),
            "tumor_state_total":         len(tumor_regulators),
            "conserved":                 reg_conserved,
            "conserved_count":           len(reg_conserved),
            "conserved_fraction":        reg_conserved_fraction,
            "population_averaged_only":  reg_pop_only,
            "tumor_state_only":          reg_tumor_only,
        },
        "targets": {
            "population_averaged_total": len(pop_targets),
            "tumor_state_total":         len(tumor_targets),
            "conserved_count":           len(tgt_conserved),
            "conserved_fraction":        tgt_conserved_fraction,
            "population_averaged_only":  tgt_pop_only,
            "tumor_state_only":          tgt_tumor_only,
        },
        "interpretation": {
            "regulatory_rewiring":              rewiring,
            "conserved_fraction_regulators":    reg_conserved_fraction,
            "tumor_specific_regulator_count":   len(reg_tumor_only),
        },
    }
