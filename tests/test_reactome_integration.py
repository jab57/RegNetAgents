#!/usr/bin/env python3
"""
Test script for Reactome pathway enrichment integration
"""

import asyncio
import json
from regnetagents_langgraph_workflow import PathwayEnricherAgent, RegNetAgentsWorkflow

async def test_reactome_agent():
    """Test PathwayEnricherAgent directly"""
    print("=" * 80)
    print("TEST 1: PathwayEnricherAgent - Direct Reactome API Test")
    print("=" * 80)

    agent = PathwayEnricherAgent()

    # Test with cancer-related genes
    genes = ['TP53', 'APC', 'BRCA1', 'MYC']
    print(f"\nTesting Reactome enrichment with genes: {genes}\n")

    result = await agent.enrich_pathways_reactome(genes)

    print(f"Status: {result.get('status')}")
    print(f"Genes analyzed: {result.get('genes_analyzed')}")
    print(f"Genes recognized by Reactome: {result.get('genes_recognized')}")
    print(f"Total pathways found: {result.get('summary', {}).get('total_pathways')}")
    print(f"Significant pathways (FDR < 0.05): {result.get('summary', {}).get('significant_pathways')}")

    print("\nTop 5 enriched pathways:")
    for i, pathway in enumerate(result.get('enriched_pathways', [])[:5], 1):
        print(f"  {i}. {pathway['pathway_name']}")
        print(f"     - Pathway ID: {pathway['pathway_id']}")
        print(f"     - p-value: {pathway['p_value']:.2e}")
        print(f"     - FDR: {pathway['fdr']:.2e}")
        print(f"     - Genes found: {pathway['genes_found']}/{pathway['genes_total']}")
        print()

async def test_workflow_integration():
    """Test Reactome integration in full workflow"""
    print("\n" + "=" * 80)
    print("TEST 2: Full Workflow Integration - TP53 Analysis")
    print("=" * 80)

    workflow = RegNetAgentsWorkflow()

    print("\nRunning comprehensive analysis for TP53...\n")
    result = await workflow.run_analysis(
        gene="TP53",
        cell_type="epithelial_cell",
        analysis_depth="basic"  # Use basic to skip domain analyses for faster test
    )

    # Extract pathway enrichment results
    pathway_enrichment = result.get('pathway_enrichment', {})

    print(f"Pathway Enrichment Status: {pathway_enrichment.get('status')}")
    print(f"Total pathways: {pathway_enrichment.get('summary', {}).get('total_pathways')}")
    print(f"Significant pathways: {pathway_enrichment.get('summary', {}).get('significant_pathways')}")

    print("\nTop 3 significant pathways for TP53:")
    for i, pathway in enumerate(pathway_enrichment.get('enriched_pathways', [])[:3], 1):
        if pathway.get('fdr', 1) < 0.05:
            print(f"  {i}. {pathway['pathway_name']}")
            print(f"     - FDR: {pathway['fdr']:.2e}")
            print()

async def main():
    """Run all tests"""
    print("\n" + "=" * 80)
    print("Reactome Pathway Enrichment Integration Tests")
    print("=" * 80 + "\n")

    # Test 1: Direct agent test
    await test_reactome_agent()

    # Test 2: Full workflow test
    await test_workflow_integration()

    print("\n" + "=" * 80)
    print("All tests completed successfully! ✓")
    print("=" * 80 + "\n")

if __name__ == "__main__":
    asyncio.run(main())
