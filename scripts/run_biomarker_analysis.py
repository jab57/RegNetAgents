#!/usr/bin/env python3
"""
Enhanced biomarker discovery analysis with full LLM insights
Saves both summary metrics AND detailed LLM-generated rationales
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

def print_sample_rationales(gene, llm_insights):
    """Print sample LLM rationales to show they're working"""
    print(f"\n  === LLM-Generated Insights for {gene} ===")

    # Clinical insights
    clinical = llm_insights.get('clinical_analysis', {})
    if clinical.get('llm_powered') and clinical.get('insights'):
        insights = clinical['insights']
        if 'biomarker_rationale' in insights:
            print(f"\n  Biomarker Rationale:")
            rationale = insights['biomarker_rationale']
            # Print first 200 chars
            print(f"    {rationale[:200]}...")
        if 'disease_rationale' in insights:
            print(f"\n  Disease Association:")
            rationale = insights['disease_rationale']
            print(f"    {rationale[:200]}...")

    # Cancer insights
    cancer = llm_insights.get('cancer_analysis', {})
    if cancer.get('llm_powered'):
        llm_rationale = cancer.get('llm_rationale', '')
        if llm_rationale and isinstance(llm_rationale, str):
            print(f"\n  Cancer Analysis Summary:")
            print(f"    {llm_rationale[:200]}...")

    # Systems insights
    systems = llm_insights.get('systems_analysis', {})
    if systems.get('llm_powered') and systems.get('insights'):
        insights = systems['insights']
        if 'centrality_rationale' in insights:
            print(f"\n  Network Centrality:")
            rationale = insights['centrality_rationale']
            print(f"    {rationale[:200]}...")

async def main():
    print("=" * 80)
    print("RegNetAgents Biomarker Discovery with Full LLM Insights")
    print("=" * 80)

    genes = ['MYC', 'CTNNB1', 'CCND1', 'TP53', 'KRAS']

    print(f"\nAnalyzing {len(genes)} candidate genes...")
    print("Genes:", ", ".join(genes))
    print("\nThis will generate:")
    print("  1. Summary metrics (biomarker_results.json)")
    print("  2. Detailed LLM insights (biomarker_llm_insights.json)")
    print("  3. Individual gene reports with rationales")

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
    detailed_insights = {}

    for gene, result in zip(genes, results):
        if isinstance(result, Exception):
            print(f"\n{gene}: FAILED - {str(result)}")
            continue

        network = result.get("network_analysis", {})
        domain_analysis = result.get("domain_analysis", {})
        cancer_info = domain_analysis.get("cancer_analysis", {})
        clinical_info = domain_analysis.get("clinical_analysis", {})
        drug_info = domain_analysis.get("drug_analysis", {})
        systems_info = domain_analysis.get("systems_analysis", {})
        pathway_info = result.get("pathway_enrichment", {})
        perturbation_data = result.get("perturbation_analysis", {})

        print(f"\n{gene}:")
        print(f"  Regulatory Role: {network.get('regulatory_role', 'unknown')}")
        print(f"  Targets: {network.get('num_targets', 0)}")
        print(f"  Regulators: {network.get('num_regulators', 0)}")

        cancer_insights = cancer_info.get("insights", {}) if cancer_info else {}
        clinical_insights = clinical_info.get("insights", {}) if clinical_info else {}
        drug_insights = drug_info.get("insights", {}) if drug_info else {}

        print(f"  Biomarker Type: {clinical_insights.get('biomarker_utility', 'unknown')}")
        print(f"  Oncogenic Potential: {cancer_insights.get('oncogenic_potential', 'unknown')}")
        print(f"  Therapeutic Score: {cancer_insights.get('therapeutic_target_score', 0)}")
        print(f"  Druggability Score: {drug_insights.get('druggability_score', 0)}")
        print(f"  Clinical Actionability: {clinical_insights.get('clinical_actionability', 'unknown')}")
        print(f"  Pathways Found: {pathway_info.get('summary', {}).get('total_pathways', 0)}")

        # Show LLM powered status
        llm_status = []
        for domain, data in domain_analysis.items():
            if data and isinstance(data, dict) and data.get('llm_powered'):
                llm_status.append(domain.replace('_analysis', ''))
        print(f"  LLM-Powered Domains: {', '.join(llm_status) if llm_status else 'None'}")

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

        # Extract LLM insights
        llm_insights = extract_llm_insights(domain_analysis)

        # Print sample rationales
        print_sample_rationales(gene, llm_insights)

        # Store detailed insights
        detailed_insights[gene] = {
            'gene': gene,
            'cell_type': 'epithelial_cell',
            'network_analysis': network,
            'llm_insights': llm_insights,
            'pathway_summary': pathway_info.get('summary', {}),
            'cross_cell_analysis': result.get('cross_cell_comparison', {}),
            'key_insights': result.get('key_insights', {}),
            'workflow_metadata': result.get('workflow_metadata', {})
        }

        # Create summary entry
        biomarker = {
            "gene": gene,
            "regulatory_role": network.get('regulatory_role', 'unknown'),
            "num_targets": network.get('num_targets', 0),
            "num_regulators": network.get('num_regulators', 0),
            "biomarker_type": clinical_insights.get('biomarker_utility', 'unknown'),
            "oncogenic_potential": cancer_insights.get('oncogenic_potential', 'unknown'),
            "therapeutic_score": cancer_insights.get('therapeutic_target_score', 0),
            "druggability_score": drug_insights.get('druggability_score', 0),
            "clinical_actionability": clinical_insights.get('clinical_actionability', 'unknown'),
            "num_pathways": pathway_info.get('summary', {}).get('total_pathways', 0),
            "llm_powered": llm_status
        }
        biomarker_summary.append(biomarker)

    # Create results directory if needed
    os.makedirs("results", exist_ok=True)

    # Save biomarker summary
    summary_output = {
        "execution_time_seconds": execution_time,
        "genes_analyzed": genes,
        "biomarker_results": biomarker_summary
    }

    summary_file = os.path.join("results", "biomarker_results.json")
    with open(summary_file, "w", encoding='utf-8') as f:
        json.dump(summary_output, f, indent=2)

    print("\n" + "=" * 80)
    print("SAVING DETAILED LLM INSIGHTS")
    print("=" * 80)

    # Save detailed insights with all LLM rationales
    insights_file = os.path.join("results", "biomarker_llm_insights.json")
    with open(insights_file, "w", encoding='utf-8') as f:
        json.dump(detailed_insights, f, indent=2)

    print(f"  [OK] biomarker_llm_insights.json ({len(detailed_insights)} genes)")

    # Save individual gene reports
    gene_reports_saved = 0
    for gene, insights in detailed_insights.items():
        filename = f"{gene.lower()}_detailed_report.json"
        filepath = os.path.join("results", filename)

        with open(filepath, "w", encoding='utf-8') as f:
            json.dump(insights, f, indent=2)

        gene_reports_saved += 1

    print(f"  [OK] {gene_reports_saved} individual gene reports")

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
    print(f"\nFiles Created:")
    print(f"  1. {summary_file} - Summary metrics")
    print(f"  2. {insights_file} - Detailed LLM insights")
    print(f"  3. {gene_reports_saved} individual gene reports")
    print(f"  4. {len(perturbation_files)} perturbation analyses")
    print("=" * 80)

if __name__ == "__main__":
    asyncio.run(main())
