#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Demo: Biomarker Panel Analysis
==============================

This script demonstrates how to use RegNetAgents programmatically WITHOUT
requiring Claude Desktop or MCP server setup. It's a great starting point
for integrating RegNetAgents into your own analysis pipelines.

Example: Colorectal Cancer Biomarker Panel
------------------------------------------
Analyzes 5 well-known colorectal cancer-associated genes:
- MYC, CTNNB1, CCND1, TP53, KRAS

What this demo shows:
---------------------
1. How to import and use the RegNetAgentsWorkflow class
2. Running parallel multi-gene analysis
3. Extracting network topology metrics
4. Therapeutic target prioritization
5. Pathway enrichment via Reactome API

Requirements:
-------------
- Python 3.8+
- pip install -r requirements.txt
- Network cache files in models/networks/ (included in repository)
- Internet connection (for Reactome pathway enrichment only)
- NO Claude Desktop required
- NO Ollama/LLM required (uses rule-based analysis by default)

Usage:
------
python demo_biomarker_analysis.py

Customize:
----------
Edit DEMO_GENES and CELL_TYPE below to analyze your own gene panel.

Expected Runtime: 15-30 seconds
"""

import asyncio
import json
import os
import sys
import time
from datetime import datetime

# Fix Windows console encoding
if sys.platform == 'win32':
    import codecs
    sys.stdout = codecs.getwriter('utf-8')(sys.stdout.buffer, 'strict')
    sys.stderr = codecs.getwriter('utf-8')(sys.stderr.buffer, 'strict')

# Import the core workflow (no MCP/Claude Desktop required)
from regnetagents_langgraph_workflow import RegNetAgentsWorkflow

# ============================================================================
# CUSTOMIZE YOUR ANALYSIS HERE
# ============================================================================

# Genes to analyze - modify this list for your own analysis
DEMO_GENES = ['MYC', 'CTNNB1', 'CCND1', 'TP53', 'KRAS']

# Cell type - choose from available networks:
# epithelial_cell, cd14_monocytes, cd16_monocytes, cd20_b_cells,
# cd4_t_cells, cd8_t_cells, erythrocytes, nk_cells, nkt_cells,
# monocyte-derived_dendritic_cells
CELL_TYPE = "epithelial_cell"

# Enable LLM-powered insights (requires Ollama with llama3.1:8b)
USE_LLM = False

# ============================================================================


def print_banner(text):
    """Print formatted banner"""
    print("\n" + "=" * 80)
    print(f"  {text}")
    print("=" * 80)


def check_ollama_available():
    """Check if Ollama is running and has the required model"""
    try:
        import ollama
        models = ollama.list()
        model_names = [m['name'] for m in models.get('models', [])]
        has_model = any('llama3.1' in name for name in model_names)
        if has_model:
            return True, "Ollama available with llama3.1"
        else:
            return False, "Ollama running but llama3.1:8b not found"
    except Exception as e:
        return False, f"Ollama not available: {e}"


def verify_environment():
    """Verify required files and environment"""
    print_banner("ENVIRONMENT CHECK")

    issues = []

    # Check network cache files
    network_dir = f"models/networks/{CELL_TYPE}"
    cache_file = os.path.join(network_dir, "network_index.pkl")

    if not os.path.exists(cache_file):
        issues.append(f"Missing network cache: {cache_file}")
    else:
        file_size = os.path.getsize(cache_file) / (1024 * 1024)  # MB
        print(f"[OK] Network cache found: {cache_file} ({file_size:.1f} MB)")

    # Check workflow module
    try:
        from regnetagents_langgraph_workflow import RegNetAgentsWorkflow
        print("[OK] Workflow module loaded")
    except ImportError as e:
        issues.append(f"Cannot import workflow: {e}")

    # Check dependencies
    try:
        import networkx
        import requests
        print("[OK] Required dependencies available")
    except ImportError as e:
        issues.append(f"Missing dependency: {e}")

    # Check Ollama (optional)
    global USE_LLM
    if USE_LLM:
        ollama_available, ollama_msg = check_ollama_available()
        if ollama_available:
            print(f"[OK] {ollama_msg}")
        else:
            print(f"[WARN] {ollama_msg}")
            print("  -> Using rule-based analysis instead")
            USE_LLM = False

    if issues:
        print("\n[FAIL] ENVIRONMENT ISSUES:")
        for issue in issues:
            print(f"  - {issue}")
        print("\nPlease run: pip install -r requirements.txt")
        sys.exit(1)

    print("\n[OK] Environment ready")


async def run_analysis():
    """Run multi-gene analysis"""
    print_banner(f"ANALYZING {len(DEMO_GENES)} GENES")

    print(f"\nGenes: {', '.join(DEMO_GENES)}")
    print(f"Cell Type: {CELL_TYPE}")
    print(f"Analysis Mode: {'LLM-powered' if USE_LLM else 'Rule-based (fast)'}")

    workflow = RegNetAgentsWorkflow()

    print("\nRunning parallel analysis...")
    start_time = time.time()

    # Run all genes in parallel
    tasks = [
        workflow.run_analysis(
            gene=gene,
            cell_type=CELL_TYPE,
            analysis_depth="comprehensive"
        )
        for gene in DEMO_GENES
    ]

    results = await asyncio.gather(*tasks, return_exceptions=True)
    execution_time = time.time() - start_time

    print(f"[OK] Analysis complete in {execution_time:.2f} seconds")

    return results, execution_time


def extract_key_metrics(gene, result):
    """Extract key metrics from analysis result"""
    if isinstance(result, Exception):
        return {
            'gene': gene,
            'status': 'FAILED',
            'error': str(result)
        }

    network = result.get('network_analysis', {})
    pathway = result.get('pathway_enrichment', {})
    targets = result.get('therapeutic_target_prioritization', {})
    domain = result.get('domain_analysis', {})

    # Extract domain insights
    cancer = domain.get('cancer_analysis', {})
    cancer_insights = cancer.get('insights', {}) if cancer else {}
    clinical = domain.get('clinical_analysis', {})
    clinical_insights = clinical.get('insights', {}) if clinical else {}
    drug = domain.get('drug_analysis', {})
    drug_insights = drug.get('insights', {}) if drug else {}

    metrics = {
        'gene': gene,
        'status': 'SUCCESS',
        'network': {
            'regulatory_role': network.get('regulatory_role', 'unknown'),
            'num_regulators': network.get('num_regulators', 0),
            'num_targets': network.get('num_targets', 0),
            'in_degree': network.get('in_degree', 0),
            'out_degree': network.get('out_degree', 0),
        },
        'pathways': {
            'total_pathways': pathway.get('summary', {}).get('total_pathways', 0),
            'significant_pathways': pathway.get('summary', {}).get('significant_pathways', 0),
        },
        'domain_insights': {
            'biomarker_type': clinical_insights.get('biomarker_utility', 'unknown'),
            'oncogenic_potential': cancer_insights.get('oncogenic_potential', 'unknown'),
            'therapeutic_score': cancer_insights.get('therapeutic_target_score', 0),
            'druggability_score': drug_insights.get('druggability_score', 0),
            'clinical_actionability': clinical_insights.get('clinical_actionability', 'unknown'),
        },
        'therapeutic_targets': None
    }

    # Add therapeutic target prioritization if available
    if targets and targets.get('ranked_regulators'):
        ranked = targets.get('ranked_regulators', [])
        top_5 = ranked[:5] if len(ranked) >= 5 else ranked

        metrics['therapeutic_targets'] = {
            'total_regulators_analyzed': len(ranked),
            'top_5_regulators': [
                {
                    'regulator': r['regulator'],
                    'pagerank': r.get('centrality_metrics', {}).get('pagerank', 0),
                    'degree_centrality': r.get('centrality_metrics', {}).get('degree_centrality', 0),
                    'regulatory_loss_pct': r.get('regulatory_loss_pct', 0),
                }
                for r in top_5
            ]
        }

    return metrics


def display_results(results_summary):
    """Display results in a readable format"""
    print_banner("ANALYSIS RESULTS")

    for result in results_summary:
        gene = result['gene']

        if result['status'] == 'FAILED':
            print(f"\n[FAIL] {gene}: {result.get('error', 'unknown error')}")
            continue

        print(f"\n[OK] {gene}:")
        print(f"  Regulatory Role: {result['network']['regulatory_role']}")
        print(f"  Regulators: {result['network']['num_regulators']}")
        print(f"  Targets: {result['network']['num_targets']}")
        print(f"  Pathways: {result['pathways']['total_pathways']}")

        if result['therapeutic_targets']:
            print(f"  Therapeutic Targets Analyzed: {result['therapeutic_targets']['total_regulators_analyzed']}")
            if result['therapeutic_targets']['top_5_regulators']:
                top_reg = result['therapeutic_targets']['top_5_regulators'][0]
                print(f"  Top Regulator: {top_reg['regulator']} (PageRank: {top_reg['pagerank']:.4f})")


def save_results(results_summary, execution_time):
    """Save results to JSON file"""
    os.makedirs("results", exist_ok=True)

    output = {
        'timestamp': datetime.now().isoformat(),
        'execution_time_seconds': execution_time,
        'genes_analyzed': DEMO_GENES,
        'cell_type': CELL_TYPE,
        'analysis_mode': 'LLM-powered' if USE_LLM else 'rule-based',
        'results': results_summary,
    }

    output_file = "results/demo_biomarker_results.json"
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(output, f, indent=2)

    print(f"\n[OK] Results saved: {output_file}")
    return output_file


async def main():
    """Main execution function"""
    print("\n" + "=" * 80)
    print("  RegNetAgents - Biomarker Panel Analysis Demo")
    print("  Standalone Example (No Claude Desktop Required)")
    print("=" * 80)

    # Step 1: Verify environment
    verify_environment()

    # Step 2: Run analysis
    results, execution_time = await run_analysis()

    # Step 3: Extract metrics
    print_banner("EXTRACTING METRICS")
    results_summary = []

    for gene, result in zip(DEMO_GENES, results):
        metrics = extract_key_metrics(gene, result)
        results_summary.append(metrics)

    print(f"[OK] Extracted metrics for {len(results_summary)} genes")

    # Step 4: Display results
    display_results(results_summary)

    # Step 5: Save results
    print_banner("SAVING RESULTS")
    output_file = save_results(results_summary, execution_time)

    # Summary
    print_banner("DEMO COMPLETE")
    print(f"\nAnalyzed {len(DEMO_GENES)} genes in {execution_time:.2f} seconds")
    print(f"Results saved to: {output_file}")
    print("\nTo customize this demo:")
    print("  1. Edit DEMO_GENES list at the top of the script")
    print("  2. Change CELL_TYPE to a different cell type")
    print("  3. Set USE_LLM = True for LLM-powered insights (requires Ollama)")
    print("\n" + "=" * 80)


if __name__ == "__main__":
    asyncio.run(main())
