"""
context_comparison.py — Cross-context regulatory comparison example.

Compares regulatory wiring for MYC across population-averaged epithelial
(GREmLN) and colorectal tumor-state (TCGA COAD) networks.

Requirements:
    TCGA caches must be built first:
        python scripts/build_tcga_cache.py --all

Usage:
    python examples/context_comparison.py
"""

import json
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from regnetagents_langgraph_workflow import RegNetAgentsWorkflow
from regnetagents.context_comparison import compare_network_contexts


def main():
    print("Loading RegNetAgents workflow...")
    workflow = RegNetAgentsWorkflow()
    agent = workflow.modeling_agent

    if agent.tcga_cache is None or not agent.tcga_cache.tcga_indices:
        print("\nTCGA cache not found. Build it first:")
        print("    python scripts/build_tcga_cache.py --all")
        sys.exit(1)

    print("\nComparing MYC regulatory wiring:")
    print("  Population-averaged context : epithelial_cell (GREmLN)")
    print("  Tumor-state context         : tcga_coad (colorectal)")
    print()

    result = compare_network_contexts(
        agent=agent,
        gene="MYC",
        cancer_type="coad",
        cell_type="epithelial_cell",
    )

    if result.get("error"):
        print(f"Error: {result['message']}")
        sys.exit(1)

    reg = result["regulators"]
    tgt = result["targets"]
    interp = result["interpretation"]

    print(f"Gene: {result['gene']}")
    print(f"Regulatory rewiring: {interp['regulatory_rewiring'].upper()}")
    print(f"Regulator conserved fraction: {interp['conserved_fraction_regulators']:.1%}")
    print()
    print("--- Regulators ---")
    print(f"  Population-averaged total : {reg['population_averaged_total']}")
    print(f"  Tumor-state total         : {reg['tumor_state_total']}")
    print(f"  Conserved ({reg['conserved_count']})             : {', '.join(reg['conserved'][:10]) or 'none'}")
    print(f"  Population-averaged only  : {', '.join(reg['population_averaged_only'][:10]) or 'none'}")
    print(f"  Tumor-state only ({interp['tumor_specific_regulator_count']})     : {', '.join(reg['tumor_state_only'][:10]) or 'none'}")
    print()
    print("--- Targets ---")
    print(f"  Population-averaged total : {tgt['population_averaged_total']}")
    print(f"  Tumor-state total         : {tgt['tumor_state_total']}")
    print(f"  Conserved                 : {tgt['conserved_count']} ({tgt['conserved_fraction']:.1%})")
    print()
    print("Full result (JSON):")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
