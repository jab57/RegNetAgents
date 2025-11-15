#!/usr/bin/env python3
"""
Simple biomarker discovery analysis - Windows compatible
Extracts both biomarker summary AND perturbation analysis results
"""
import asyncio
import json
import os
import sys
import time

# Add parent directory to path so we can import the workflow
script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(script_dir)
sys.path.insert(0, parent_dir)

from regnetagents_langgraph_workflow import RegNetAgentsWorkflow

async def main():
    print("=" * 80)
    print("RegNetAgents Biomarker Discovery for Colorectal Cancer")
    print("=" * 80)

    genes = ['MYC', 'CTNNB1', 'CCND1', 'TP53', 'KRAS']

    print(f"\nAnalyzing {len(genes)} candidate genes...")
    print("Genes:", ", ".join(genes))

    workflow = RegNetAgentsWorkflow()

    start_time = time.time()

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

    biomarker_summary = []
    perturbation_results = {}

    for gene, result in zip(genes, results):
        if isinstance(result, Exception):
            print(f"\n{gene}: FAILED - {str(result)}")
            continue

        network = result.get("network_analysis", {})
        cancer_info = result.get("domain_analysis", {}).get("cancer_analysis", {})
        clinical_info = result.get("domain_analysis", {}).get("clinical_analysis", {})
        pathway_info = result.get("pathway_enrichment", {})
        perturbation_data = result.get("perturbation_analysis", {})

        print(f"\n{gene}:")
        print(f"  Regulatory Role: {network.get('regulatory_role', 'unknown')}")
        print(f"  Targets: {network.get('num_targets', 0)}")
        print(f"  Regulators: {network.get('num_regulators', 0)}")

        cancer_insights = cancer_info.get("insights", {}) if cancer_info else {}
        clinical_insights = clinical_info.get("insights", {}) if clinical_info else {}

        print(f"  Biomarker Type: {clinical_insights.get('biomarker_utility', 'unknown')}")
        print(f"  Oncogenic Potential: {cancer_insights.get('oncogenic_potential', 'unknown')}")
        print(f"  Therapeutic Score: {cancer_insights.get('therapeutic_target_score', 0)}")
        print(f"  Clinical Actionability: {clinical_insights.get('clinical_actionability', 'unknown')}")
        print(f"  Pathways Found: {pathway_info.get('summary', {}).get('total_pathways', 0)}")

        # Extract perturbation analysis if available
        if perturbation_data:
            num_perturbations = len(perturbation_data.get('perturbation_results', []))
            print(f"  Perturbation Analysis: {num_perturbations} regulators analyzed")

            # Get top regulator by PageRank
            if perturbation_data.get('rankings', {}).get('by_pagerank'):
                top = perturbation_data['rankings']['by_pagerank'][0]
                print(f"  Top Regulator (PageRank): {top['regulator']} ({top['score']:.4f})")
                perturbation_results[gene] = perturbation_data
        else:
            print(f"  Perturbation Analysis: Not available (insufficient regulators)")

        biomarker = {
            "gene": gene,
            "regulatory_role": network.get('regulatory_role', 'unknown'),
            "num_targets": network.get('num_targets', 0),
            "num_regulators": network.get('num_regulators', 0),
            "biomarker_type": clinical_insights.get('biomarker_utility', 'unknown'),
            "oncogenic_potential": cancer_insights.get('oncogenic_potential', 'unknown'),
            "therapeutic_score": cancer_insights.get('therapeutic_target_score', 0),
            "clinical_actionability": clinical_insights.get('clinical_actionability', 'unknown'),
            "num_pathways": pathway_info.get('summary', {}).get('total_pathways', 0)
        }
        biomarker_summary.append(biomarker)

    # Save biomarker summary
    output = {
        "execution_time_seconds": execution_time,
        "genes_analyzed": genes,
        "biomarker_results": biomarker_summary
    }

    # Create results directory if needed
    os.makedirs("results", exist_ok=True)

    summary_file = os.path.join("results", "biomarker_results.json")
    with open(summary_file, "w") as f:
        json.dump(output, f, indent=2)

    print("\n" + "=" * 80)
    print("SAVING PERTURBATION RESULTS")
    print("=" * 80)

    # Save individual perturbation results
    perturbation_files = []
    for gene, perturb_data in perturbation_results.items():
        filename = f"{gene.lower()}_perturbation_standard_centrality.json"
        filepath = os.path.join("results", filename)

        with open(filepath, "w", encoding='utf-8') as f:
            json.dump(perturb_data, f, indent=2)

        print(f"  [OK] {filename}")
        perturbation_files.append(filename)

    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Genes Analyzed: {len(genes)}")
    print(f"Execution Time: {execution_time:.2f} seconds")
    print(f"Biomarker Summary: {summary_file}")
    print(f"Perturbation Results: {len(perturbation_files)} files saved to results/")
    print("=" * 80)

if __name__ == "__main__":
    asyncio.run(main())
