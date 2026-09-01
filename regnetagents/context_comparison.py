"""
context_comparison.py — Compare regulatory wiring across network contexts.

Computes the overlap between population-averaged (GREmLN) and tumor-state (TCGA)
regulatory programs for a given gene. All logic is rule-based — no LLM required.
"""

from __future__ import annotations
from typing import TYPE_CHECKING

from . import driver_gene_client

if TYPE_CHECKING:
    from regnetagents_langgraph_workflow import RegNetAgentsModelingAgent


def _extract_symbols_with_weights(entries) -> tuple[set, dict]:
    """Collapse gene_neighbors entries to a symbol set plus a symbol->weight map.

    Entries are dicts with a "gene" key and (for TCGA-sourced entries) a
    "likelihood" MI-weight key; legacy/non-dict entries (plain symbol strings)
    yield an empty weights dict.
    """
    if not entries:
        return set(), {}
    if isinstance(entries[0], dict):
        symbols = {e["gene"] for e in entries if "gene" in e}
        weights = {e["gene"]: e.get("likelihood") for e in entries if "gene" in e}
        return symbols, weights
    return set(entries), {}


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
    # Both GREmLN and TCGA gene_neighbors return regulators/targets at the
    # top level of the result dict (not nested under a "gene_neighbors" key).
    pop_regulators, _                    = _extract_symbols_with_weights(pop_result.get("regulators", []))
    pop_targets, _                       = _extract_symbols_with_weights(pop_result.get("targets", []))
    tumor_regulators, tumor_reg_weights  = _extract_symbols_with_weights(tumor_result.get("regulators", []))
    tumor_targets, _                     = _extract_symbols_with_weights(tumor_result.get("targets", []))

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

    # --- Cancer-driver annotation (IntOGen, always-on) ----------------------
    # In-memory dict lookup after the first call warms the cache; annotates the
    # union of all regulators. `driver_gene_roles` maps every regulator to its
    # collapsed IntOGen role or None (non-driver). `tumor_state_only_known_drivers`
    # is the "signal vs. noise" list: tumor-acquired regulators that are known
    # cancer drivers -- filter on the *value* (a non-None role), never on `g in
    # driver_gene_roles` (every regulator is a key of that dict).
    driver_roles_map = driver_gene_client.get_driver_roles()
    driver_gene_roles = {
        g: driver_roles_map.get(g)
        for g in (reg_conserved + reg_pop_only + reg_tumor_only)
    }
    tumor_state_only_known_drivers = sorted(
        g for g in reg_tumor_only if driver_gene_roles.get(g) is not None
    )
    # Tissue-matched subset: tumor-acquired known drivers that IntOGen called
    # specifically in *this* cancer_type (not just somewhere pan-cancer).
    # Presence only -- `driver_gene_roles` values stay pan-cancer regardless.
    tumor_state_only_tissue_matched_drivers = sorted(
        g for g in tumor_state_only_known_drivers
        if driver_gene_client.is_tissue_matched(g, cancer_type)
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
            "tumor_state_only_weights":  {g: tumor_reg_weights.get(g) for g in reg_tumor_only},
            "driver_gene_roles":              driver_gene_roles,
            "tumor_state_only_known_drivers": tumor_state_only_known_drivers,
            "tumor_state_only_tissue_matched_drivers": tumor_state_only_tissue_matched_drivers,
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
            "tumor_state_only_known_driver_count": len(tumor_state_only_known_drivers),
            "tumor_state_only_tissue_matched_driver_count": len(tumor_state_only_tissue_matched_drivers),
        },
        "driver_annotation_available": driver_gene_client.driver_data_available(),
    }
