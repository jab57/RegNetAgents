"""
Verify pathway enrichment counts for the 5-gene CRC panel.
Runs actual Reactome API calls to get real pathway counts.
"""

import json
import asyncio
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from regnetagents_langgraph_workflow import PathwayEnricherAgent, RegNetAgentsModelingAgent, RegNetAgentsCache

# Global cache (initialize once)
CACHE = None

async def get_pathway_count(gene_symbol, cell_type="epithelial_cell"):
    """Get actual pathway enrichment count for a gene."""

    global CACHE
    if CACHE is None:
        print("Initializing network cache...")
        CACHE = RegNetAgentsCache()

    # Initialize agents
    network_agent = RegNetAgentsModelingAgent(CACHE)
    pathway_enricher = PathwayEnricherAgent()

    # Get network analysis
    print(f"\n{'='*60}")
    print(f"Analyzing {gene_symbol}...")
    print(f"{'='*60}")

    from regnetagents_langgraph_workflow import CellType
    cell_type_enum = CellType.EPITHELIAL_CELLS
    network_result = await network_agent.analyze_gene_network_context(gene_symbol, cell_type_enum)

    if not network_result['in_network']:
        print(f"[ERROR] {gene_symbol} not found in network")
        return None

    num_regulators = network_result['num_regulators']
    num_targets = network_result['num_targets']

    print(f"  Regulators: {num_regulators}")
    print(f"  Targets: {num_targets}")

    # Build gene list for pathway enrichment
    # Gene + top 10 regulators + top 10 targets
    gene_list = [gene_symbol]

    # Add regulators (up to 10)
    if network_result.get('regulators'):
        regulator_symbols = [network_agent._convert_ensembl_to_symbol(ens_id)
                            for ens_id in network_result['regulators'][:10]]
        gene_list.extend([s for s in regulator_symbols if s and s != "Unknown"])

    # Add targets (up to 10)
    if network_result.get('targets'):
        target_symbols = [network_agent._convert_ensembl_to_symbol(ens_id)
                         for ens_id in network_result['targets'][:10]]
        gene_list.extend([s for s in target_symbols if s and s != "Unknown"])

    print(f"  Gene list size: {len(gene_list)} genes")
    print(f"  Gene list: {', '.join(gene_list[:15])}{'...' if len(gene_list) > 15 else ''}")

    # Run pathway enrichment
    print(f"  Querying Reactome API...")
    pathway_result = await pathway_enricher.enrich_pathways_reactome(gene_list)

    # Debug: Print API response structure
    print(f"  DEBUG: API response keys: {list(pathway_result.keys()) if pathway_result else 'None'}")
    print(f"  DEBUG: Status: {pathway_result.get('status')}")
    if pathway_result.get('status') == 'success':
        print(f"  DEBUG: Genes recognized: {pathway_result.get('genes_recognized', 0)}")
        print(f"  DEBUG: Total enriched pathways: {pathway_result.get('summary', {}).get('total_pathways', 0)}")
        if pathway_result.get('enriched_pathways'):
            print(f"  DEBUG: First pathway: {pathway_result['enriched_pathways'][0] if pathway_result['enriched_pathways'] else 'None'}")

    # Count significant pathways (FDR < 0.05)
    num_pathways = 0
    if pathway_result and pathway_result.get('status') == 'success':
        num_pathways = len(pathway_result.get('enriched_pathways', []))

    print(f"  [OK] Significant pathways (FDR < 0.05): {num_pathways}")

    return {
        'gene': gene_symbol,
        'num_regulators': num_regulators,
        'num_targets': num_targets,
        'gene_list_size': len(gene_list),
        'num_pathways': num_pathways
    }

async def main():
    """Run pathway analysis for all 5 genes."""

    genes = ['MYC', 'CTNNB1', 'CCND1', 'TP53', 'KRAS']

    print("\n" + "="*60)
    print("PATHWAY ENRICHMENT VERIFICATION")
    print("Running Reactome API queries for 5-gene CRC panel")
    print("="*60)

    results = []
    for gene in genes:
        result = await get_pathway_count(gene)
        if result:
            results.append(result)
        await asyncio.sleep(1)  # Rate limiting

    # Summary table
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"{'Gene':<10} {'Regs':<6} {'Tgts':<6} {'List':<6} {'Pathways':<10}")
    print("-" * 60)
    for r in results:
        print(f"{r['gene']:<10} {r['num_regulators']:<6} {r['num_targets']:<6} "
              f"{r['gene_list_size']:<6} {r['num_pathways']:<10}")

    # Save results
    output_file = 'results/pathway_verification.json'
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)

    print(f"\n[OK] Results saved to {output_file}")

    return results

if __name__ == "__main__":
    results = asyncio.run(main())
