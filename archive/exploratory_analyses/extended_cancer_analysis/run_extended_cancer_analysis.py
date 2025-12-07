#!/usr/bin/env python3
"""
Extended multi-cancer gene analysis for manuscript validation
Analyzes genes from breast, lung, prostate, and pancreatic cancers
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

def extract_llm_insights(domain_analysis):
    """Extract all LLM rationales and insights from domain analysis"""
    llm_insights = {}

    for domain, data in domain_analysis.items():
        if not data or not isinstance(data, dict):
            continue

        domain_insights = {
            'llm_powered': data.get('llm_powered', False),
            'summary': data.get('summary', ''),
            'insights': data.get('insights', {}),
            'llm_rationale': data.get('llm_rationale', '')
        }

        # Extract domain-specific additional data
        if domain == 'cancer_analysis':
            domain_insights['cancer_pathways'] = data.get('cancer_pathways', [])
            domain_insights['biomarker_potential'] = data.get('biomarker_potential', {})
        elif domain == 'drug_analysis':
            domain_insights['cascade_effects'] = data.get('cascade_effects', {})
            domain_insights['clinical_trial_readiness'] = data.get('clinical_trial_readiness', {})
        elif domain == 'clinical_analysis':
            domain_insights['tissue_specificity'] = data.get('tissue_specificity', [])
            domain_insights['diagnostic_potential'] = data.get('diagnostic_potential', {})
        elif domain == 'systems_analysis':
            domain_insights['network_effects'] = data.get('network_effects', {})
            domain_insights['perturbation_impact'] = data.get('perturbation_impact', {})

        llm_insights[domain] = domain_insights

    return llm_insights

async def main():
    print("=" * 80)
    print("RegNetAgents Extended Multi-Cancer Validation")
    print("=" * 80)

    # Extended gene panel across multiple cancer types
    genes = {
        'Breast Cancer': ['BRCA1', 'BRCA2', 'ERBB2', 'ESR1'],
        'Lung Cancer': ['EGFR', 'ALK', 'RET'],
        'Prostate Cancer': ['AR', 'PTEN'],
        'Pancreatic Cancer': ['CDKN2A']
    }

    # Flatten to single list for analysis
    all_genes = []
    for cancer_type, gene_list in genes.items():
        all_genes.extend(gene_list)

    print(f"\nAnalyzing {len(all_genes)} genes across 4 cancer types...")
    print("\nGene Panel:")
    for cancer_type, gene_list in genes.items():
        print(f"  {cancer_type}: {', '.join(gene_list)}")

    print("\nThis will generate:")
    print("  1. Individual gene reports with network analysis")
    print("  2. Therapeutic target prioritization (where applicable)")
    print("  3. Summary table for manuscript")

    workflow = RegNetAgentsWorkflow()

    start_time = time.time()

    tasks = [
        workflow.run_analysis(
            gene=gene,
            cell_type="epithelial_cell",
            analysis_depth="comprehensive"
        )
        for gene in all_genes
    ]

    results = await asyncio.gather(*tasks, return_exceptions=True)

    execution_time = time.time() - start_time

    print(f"\nAnalysis complete in {execution_time:.2f} seconds")
    print("\n" + "=" * 80)
    print("RESULTS")
    print("=" * 80)

    validation_summary = []
    detailed_insights = {}
    perturbation_results = {}

    for gene, result in zip(all_genes, results):
        if isinstance(result, Exception):
            print(f"\n{gene}: FAILED - {str(result)}")
            validation_summary.append({
                'gene': gene,
                'status': 'failed',
                'error': str(result)
            })
            continue

        # Handle None results
        if result is None:
            print(f"\n{gene}: FAILED - No result returned")
            validation_summary.append({
                'gene': gene,
                'status': 'failed',
                'error': 'No result returned'
            })
            continue

        network = result.get("network_analysis", {})
        domain_analysis = result.get("domain_analysis", {})
        perturbation_data = result.get("perturbation_analysis", {})

        print(f"\n{gene}:")
        print(f"  Regulatory Role: {network.get('regulatory_role', 'unknown')}")
        print(f"  Targets: {network.get('num_targets', 0)}")
        print(f"  Regulators: {network.get('num_regulators', 0)}")

        # Extract perturbation analysis
        top_regulator = "N/A"
        top_pagerank = 0.0
        if perturbation_data:
            if perturbation_data.get('rankings', {}).get('by_pagerank'):
                top = perturbation_data['rankings']['by_pagerank'][0]
                top_regulator = top['regulator']
                top_pagerank = top['score']
                print(f"  Top Regulator: {top_regulator} (PageRank: {top_pagerank:.3f})")
                perturbation_results[gene] = perturbation_data

        # Extract LLM insights
        llm_insights = extract_llm_insights(domain_analysis)

        # Store detailed insights
        pathway_data = result.get('pathway_enrichment') or {}
        detailed_insights[gene] = {
            'gene': gene,
            'cell_type': 'epithelial_cell',
            'network_analysis': network,
            'llm_insights': llm_insights,
            'pathway_summary': pathway_data.get('summary', {}),
            'cross_cell_analysis': result.get('cross_cell_comparison', {}),
            'key_insights': result.get('key_insights', {}),
            'workflow_metadata': result.get('workflow_metadata', {})
        }

        # Create validation summary entry
        validation_summary.append({
            'gene': gene,
            'regulatory_role': network.get('regulatory_role', 'unknown'),
            'num_targets': network.get('num_targets', 0),
            'num_regulators': network.get('num_regulators', 0),
            'top_regulator': top_regulator,
            'top_pagerank': top_pagerank,
            'status': 'success'
        })

    # Create results directory if needed
    os.makedirs("results", exist_ok=True)

    print("\n" + "=" * 80)
    print("SAVING RESULTS")
    print("=" * 80)

    # Save validation summary
    summary_file = os.path.join("results", "extended_cancer_validation.json")
    with open(summary_file, "w", encoding='utf-8') as f:
        json.dump({
            'execution_time_seconds': execution_time,
            'genes_analyzed': all_genes,
            'cancer_types': genes,
            'validation_results': validation_summary
        }, f, indent=2)
    print(f"  [OK] extended_cancer_validation.json")

    # Save individual gene reports
    gene_reports_saved = 0
    for gene, insights in detailed_insights.items():
        filename = f"{gene.lower()}_detailed_report.json"
        filepath = os.path.join("results", filename)

        with open(filepath, "w", encoding='utf-8') as f:
            json.dump(insights, f, indent=2)

        gene_reports_saved += 1

    print(f"  [OK] {gene_reports_saved} individual gene reports")

    # Save individual perturbation results
    perturbation_files = 0
    for gene, perturb_data in perturbation_results.items():
        filename = f"{gene.lower()}_perturbation_standard_centrality.json"
        filepath = os.path.join("results", filename)

        with open(filepath, "w", encoding='utf-8') as f:
            json.dump(perturb_data, f, indent=2)

        perturbation_files += 1

    print(f"  [OK] {perturbation_files} perturbation analyses")

    print("\n" + "=" * 80)
    print("MANUSCRIPT TABLE PREVIEW")
    print("=" * 80)
    print("\nTable: Extended Multi-Cancer Gene Validation\n")
    print(f"{'Gene':<8} {'Cancer Type':<18} {'Role':<18} {'Targets':<8} {'Regs':<6} {'Top Regulator':<15}")
    print("-" * 90)

    for cancer_type, gene_list in genes.items():
        for gene in gene_list:
            summary = next((s for s in validation_summary if s['gene'] == gene), None)
            if summary and summary['status'] == 'success':
                print(f"{gene:<8} {cancer_type:<18} {summary['regulatory_role']:<18} "
                      f"{summary['num_targets']:<8} {summary['num_regulators']:<6} {summary['top_regulator']:<15}")

    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Genes Analyzed: {len(all_genes)}")
    print(f"Successful: {sum(1 for s in validation_summary if s['status'] == 'success')}")
    print(f"Failed: {sum(1 for s in validation_summary if s['status'] == 'failed')}")
    print(f"Execution Time: {execution_time:.2f} seconds")
    print(f"\nNext Steps:")
    print(f"  1. Review results in: {summary_file}")
    print(f"  2. Check individual gene reports in results/ folder")
    print(f"  3. Validate top regulators against literature")
    print("=" * 80)

if __name__ == "__main__":
    asyncio.run(main())
