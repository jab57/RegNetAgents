#!/usr/bin/env python3
"""
Generate perturbation analysis results for all 5 biomarker genes
Uses the existing RegNetAgents workflow to analyze MYC, CTNNB1, CCND1, TP53, KRAS
and save their perturbation results to JSON files.
"""

import asyncio
import json
import os
from regnetagents_langgraph_workflow import RegNetAgentsWorkflow


async def analyze_gene_perturbations(gene, cell_type='epithelial_cell'):
    """Run comprehensive analysis for a gene and extract perturbation results"""
    print(f"\n{'='*60}")
    print(f"Analyzing {gene}")
    print(f"{'='*60}")

    workflow = RegNetAgentsWorkflow()

    try:
        result = await workflow.run_analysis(
            gene=gene,
            cell_type=cell_type,
            analysis_depth='focused'  # Focused includes perturbation without Reactome delay
        )

        if result.get('status') == 'error':
            print(f"  ❌ Error: {result.get('error', 'Unknown error')}")
            return None

        # Extract perturbation results
        perturbation_results = result.get('perturbation_analysis')

        if perturbation_results:
            num_regulators = len(perturbation_results.get('perturbation_results', []))
            print(f"  ✓ Found {num_regulators} regulators")

            # Get top regulator by PageRank
            if perturbation_results.get('rankings', {}).get('by_pagerank'):
                top = perturbation_results['rankings']['by_pagerank'][0]
                print(f"  ✓ Top by PageRank: {top['regulator']} ({top['score']:.4f})")

            return perturbation_results
        else:
            print(f"  ⚠ No perturbation results (may have <5 regulators)")
            return None

    except Exception as e:
        print(f"  ❌ Exception: {type(e).__name__}: {e}")
        return None


async def main():
    """Generate perturbation results for all 5 biomarker genes"""
    print("="*60)
    print("RegNetAgents Biomarker Perturbation Analysis Generator")
    print("="*60)
    print("\nGenerating perturbation analysis for CRC biomarker panel:")
    print("  MYC, CTNNB1, CCND1, TP53, KRAS")

    genes = ['MYC', 'CTNNB1', 'CCND1', 'TP53', 'KRAS']
    results = {}

    for gene in genes:
        perturbation_data = await analyze_gene_perturbations(gene)
        if perturbation_data:
            results[gene] = perturbation_data

    # Save individual results
    print(f"\n{'='*60}")
    print("Saving Results")
    print(f"{'='*60}")

    results_dir = 'results'
    os.makedirs(results_dir, exist_ok=True)

    for gene, data in results.items():
        filename = f"{gene.lower()}_perturbation_standard_centrality.json"
        filepath = os.path.join(results_dir, filename)

        with open(filepath, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=2)

        print(f"  ✓ Saved {filepath}")

    # Create summary table
    print(f"\n{'='*60}")
    print("Summary Table")
    print(f"{'='*60}")
    print(f"{'Gene':<10} {'Regulators':<12} {'Top by PageRank':<20} {'Score':<10}")
    print("-" * 60)

    for gene in genes:
        if gene in results:
            data = results[gene]
            num_reg = len(data.get('perturbation_results', []))

            if data.get('rankings', {}).get('by_pagerank'):
                top = data['rankings']['by_pagerank'][0]
                top_reg = top['regulator']
                score = top['score']
                print(f"{gene:<10} {num_reg:<12} {top_reg:<20} {score:<10.4f}")
            else:
                print(f"{gene:<10} {num_reg:<12} {'N/A':<20} {'N/A':<10}")
        else:
            print(f"{gene:<10} {'ERROR':<12} {'N/A':<20} {'N/A':<10}")

    print(f"\n{'='*60}")
    print(f"SUCCESS! Generated {len(results)}/5 perturbation analyses")
    print(f"{'='*60}")
    print(f"\nResults saved to: {results_dir}/")
    print("Next steps:")
    print("  1. Review the perturbation results")
    print("  2. Update publication tables with top PageRank regulators")
    print("  3. Regenerate figures and poster")


if __name__ == "__main__":
    asyncio.run(main())
