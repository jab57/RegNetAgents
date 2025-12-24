#!/usr/bin/env python3
"""
Analyze cervical cancer genes for supplementary validation
Focuses on 7 genes with confirmed regulatory network coverage
"""
import asyncio
import json
import os
import sys
import time

# Add parent directory to path
script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(script_dir)
sys.path.insert(0, parent_dir)

from regnetagents_langgraph_workflow import RegNetAgentsWorkflow

async def main():
    print("=" * 80)
    print("Cervical Cancer Gene Analysis - Supplementary Validation")
    print("=" * 80)

    # 7 genes with confirmed good network coverage
    genes = ['TP53', 'BRCA1', 'BRCA2', 'MYC', 'KRAS', 'EGFR', 'CCND1']

    print(f"\nAnalyzing {len(genes)} cervical cancer genes...")
    print(f"Genes: {', '.join(genes)}")
    print("\nAll genes have confirmed regulatory network coverage (>=5 regulators)")

    workflow = RegNetAgentsWorkflow()
    start_time = time.time()

    # Run analyses in parallel
    tasks = [
        workflow.run_analysis(
            gene=gene,
            cell_type="epithelial_cell",
            analysis_depth="comprehensive"
        )
        for gene in genes
    ]

    results = await asyncio.gather(*tasks, return_exceptions=True)
    execution_time = time.time() - start_time

    print(f"\nAnalysis complete in {execution_time:.2f} seconds")
    print("\n" + "=" * 80)
    print("RESULTS")
    print("=" * 80)

    summary_data = []

    for gene, result in zip(genes, results):
        if isinstance(result, Exception):
            print(f"\n{gene}: FAILED - {str(result)}")
            summary_data.append({
                'gene': gene,
                'status': 'failed',
                'error': str(result)
            })
            continue

        network = result.get("network_analysis", {})
        target_prioritization = result.get("therapeutic_target_prioritization", {})
        pathways = result.get("pathway_enrichment", {})

        print(f"\n{gene}:")
        print(f"  Role: {network.get('regulatory_role', 'unknown')}")
        print(f"  Regulators: {network.get('num_regulators', 0)}")
        print(f"  Targets: {network.get('num_targets', 0)}")

        # Get top regulator
        top_regulator = "N/A"
        top_pagerank = 0.0
        if target_prioritization and target_prioritization.get('rankings', {}).get('by_pagerank'):
            top = target_prioritization['rankings']['by_pagerank'][0]
            top_regulator = top['regulator']
            top_pagerank = top['score']
            print(f"  Top Regulator: {top_regulator} (PageRank: {top_pagerank:.3f})")

        # Get top pathways
        top_pathways = []
        if pathways and pathways.get('significant_pathways'):
            top_pathways = [
                p['pathway_name']
                for p in pathways['significant_pathways'][:3]
            ]
            print(f"  Top Pathways: {', '.join(top_pathways[:2])}")

        # Save individual gene report
        gene_report_file = f"results/cervical_{gene.lower()}_analysis.json"
        with open(gene_report_file, 'w') as f:
            json.dump(result, f, indent=2)
        print(f"  Report saved: {gene_report_file}")

        summary_data.append({
            'gene': gene,
            'regulatory_role': network.get('regulatory_role', 'unknown'),
            'num_regulators': network.get('num_regulators', 0),
            'num_targets': network.get('num_targets', 0),
            'top_regulator': top_regulator,
            'top_pagerank': top_pagerank,
            'top_pathways': top_pathways,
            'status': 'success'
        })

    # Save summary
    summary_file = "results/cervical_cancer_supplementary_validation.json"
    summary = {
        'analysis_date': time.strftime('%Y-%m-%d'),
        'execution_time_seconds': execution_time,
        'genes_analyzed': genes,
        'cell_type': 'epithelial_cell',
        'purpose': 'Supplementary validation for cervical cancer',
        'results': summary_data
    }

    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)

    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Genes analyzed: {len(genes)}")
    print(f"Successful: {len([r for r in summary_data if r['status'] == 'success'])}")
    print(f"Failed: {len([r for r in summary_data if r['status'] == 'failed'])}")
    print(f"Execution time: {execution_time:.2f} seconds")
    print(f"\nSummary saved: {summary_file}")
    print(f"Individual reports: results/cervical_[gene]_analysis.json")
    print("=" * 80)

if __name__ == "__main__":
    asyncio.run(main())
