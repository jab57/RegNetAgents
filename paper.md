---
title: 'RegNetAgents: A Multi-Agent LLM System for Gene Regulatory Network Analysis'
tags:
  - Python
  - multi-agent systems
  - gene regulatory networks
  - LLM orchestration
  - bioinformatics
  - workflow automation
authors:
  - name: Jose A. Bird
    orcid: 0009-0006-2744-0606
    affiliation: 1
affiliations:
  - name: Bird AI Solutions
    index: 1
date: "05 February 2026"
bibliography: paper.bib
---

# Summary

RegNetAgents is a Python package for automated downstream analysis of pre-computed gene regulatory networks (GRNs) using a multi-agent architecture built on LangGraph [@langgraph]. The software coordinates four specialized domain agents (cancer biology, drug discovery, clinical relevance, systems biology) to analyze regulatory relationships from ARACNe-inferred networks [@margolin2006aracne; @lachmann2016aracne]. It can be used programmatically via Python or through MCP-compatible clients (Claude Desktop, Cursor, Zed) [@mcp] via natural language. RegNetAgents is intended for computational biologists, bioinformaticians, and researchers who need automated GRN analysis.

The architecture separates workflow logic from LLM interpretation: a fixed directed acyclic graph (DAG) controls execution order and data flow, while an optional LLM layer adds natural-language narrative rationale to deterministic rule-based domain assessments. This separation ensures that identical inputs produce identical computational outputs, with LLM variability isolated to the interpretation layer. RegNetAgents thus provides a reproducible entry point for downstream GRN interpretation.

Key features:

- **Hybrid architecture**: Combines optional LLM-powered interpretation with a deterministic rule-based layer for reproducible, auditable assessments
- **Parallel agents**: Four domain agents run concurrently for fast analysis
- **Therapeutic prioritization**: Deterministic PageRank and centrality-based ranking of upstream regulators [@mora2021effects]
- **Pathway enrichment**: Automated Reactome API queries with FDR correction [@gillespie2022reactome]

# Statement of Need

Researchers analyzing gene regulatory networks typically query multiple databases manually (STRING, BioGRID, Reactome), export data across platforms, and synthesize findings—a process requiring hours per gene and programming expertise. Established tools address individual components of this workflow: Cytoscape [@shannon2003cytoscape] provides interactive network visualization, pySCENIC [@aibar2017scenic] and CellOracle [@kamimoto2023dissecting] infer regulatory networks from expression data, and STRING [@szklarczyk2023string] catalogs known protein interactions. However, these tools operate on separate steps that must be manually connected, require programming expertise, and produce outputs that need expert interpretation across multiple biological contexts. These tools focus on network inference or visualization, whereas RegNetAgents operates downstream, automating multi-perspective interpretation of pre-computed networks. Three contrasts define RegNetAgents' novelty: (1) existing tools target network inference or visualization, while RegNetAgents targets downstream biological interpretation of pre-computed networks; (2) general-purpose agent frameworks rely on non-deterministic autonomous execution, while RegNetAgents uses a fixed DAG that guarantees identical outputs for identical inputs; and (3) individual tools each address a single analytical context, while RegNetAgents simultaneously synthesizes cancer biology, drug discovery, clinical relevance, and systems biology perspectives in one pass. **No existing tool provides an integrated, conversational workflow combining multi-domain GRN interpretation with deterministic reproducibility.**

RegNetAgents addresses this gap by combining a natural-language interface with a reproducible, multi-agent workflow that performs regulator/target identification, therapeutic prioritization, pathway enrichment, and domain-specific interpretation in a single automated pipeline, without requiring users to write code. In a five-gene colorectal cancer case study, the system analyzed 99 upstream regulators in 15–62 seconds; single-gene comprehensive analyses typically complete in under a minute on a standard laptop. The equivalent manual workflow—querying regulators and targets across 10 cell types, performing pathway enrichment for each, and synthesizing findings across four biological perspectives per gene—requires an estimated 8–16 hours, a speedup of three orders of magnitude.

The software processes pre-computed ARACNe networks from the GREmLN project [@zhang2025gremln], covering 10 cell types—including epithelial, stromal, immune, and endothelial lineages relevant to colorectal cancer biology—derived from large single-cell RNA-seq datasets via the CELLxGENE Data Portal [@megill2021cellxgene]. ARACNe is a well-validated mutual-information method that produces high-confidence transcription factor regulons by applying the data-processing inequality to large expression datasets; the GREmLN networks represent a curated, cell-type-resolved regulatory atlas validated on the colorectal tissue atlas. Unlike general-purpose agent frameworks [@autogpt; @li2023camel; @hong2023metagpt], RegNetAgents provides domain-specific scientific workflow integration with deterministic fallback guarantees for reproducible research.

# Implementation

RegNetAgents requires Python 3.10+ and uses NetworkX [@hagberg2008networkx] for graph algorithms. The workflow is implemented as a LangGraph DAG, which enforces explicit stepwise execution and ensures reproducibility across runs. This design intentionally reduces agent autonomy—agents cannot dynamically alter the workflow—in exchange for deterministic, auditable execution. Three design trade-offs follow from this choice: the DAG constraints limit agent autonomy but guarantee reproducibility and auditability; the dual-mode architecture (LLM + rule-based) trades richer interpretation for guaranteed determinism when needed; and restricting input to pre-computed ARACNe networks prioritizes validated, curated data over generality. The pipeline executes in seven deterministic steps: (1) **gene validation**—symbol normalization and network membership check; (2) **network lookup**—regulator and target retrieval per cell type; (3) **therapeutic prioritization**—PageRank- and centrality-based regulator ranking (conditional on >5 regulators); (4) **pathway enrichment**—Reactome API queries with FDR correction; (5) **cross-cell comparison**—multi-cell regulatory role and pathway consistency (conditional for hub and master regulators); (6) **parallel domain analysis**—concurrent execution of four specialized agents; and (7) **report synthesis**—cross-domain contradiction detection and narrative generation (\autoref{fig:architecture}).

![RegNetAgents multi-agent architecture. The LangGraph DAG enforces a fixed execution order: gene validation, network lookup, therapeutic prioritization, pathway enrichment, cross-cell comparison (conditional for hub/master regulators), and parallel domain analysis, followed by report synthesis. The rule-based layer produces categorical evidence assessments (high/moderate/low); the optional LLM layer adds natural-language narrative rationale without altering the assessments. This figure depicts the \texttt{comprehensive\_gene\_analysis} pipeline; additional MCP tools including \texttt{find\_master\_regulators} operate as direct network queries outside the DAG.\label{fig:architecture}](figure1_architecture.png)

Therapeutic prioritization uses NetworkX's deterministic PageRank and degree-based centrality metrics to rank upstream regulators by influence within the ARACNe-derived subnetwork. The four domain agents each produce an independent evidence assessment (high/moderate/low) for a specific biological dimension: oncogenic potential (cancer agent), druggability (drug discovery agent), clinical actionability (clinical relevance agent), and network centrality (systems biology agent). The rule-based layer—active by default (`USE_LLM_AGENTS=false`)—derives these assessments from network topology using empirically derived, cell-type-specific thresholds. For each cell-type network, the preprocessing pipeline automatically computes the 90th percentile of target count (out-degree), regulator count (in-degree), and normalized PageRank—measured across all nodes (transcription factors and their targets) in the cell-type GRN—as the **high** threshold, and the 75th percentile as the **moderate** threshold. These percentile cutoffs are stored per cell type and loaded at runtime, so the cancer agent's criterion for "high oncogenic potential" adapts to the regulatory scale of each network rather than using a fixed value. **High** is assigned when a gene satisfies two or more domain-specific criteria; **moderate** when one criterion is met; **low** when none apply. The multi-criteria requirement for **high** promotes conservative classification by requiring independent convergent signals. These thresholds are recomputed automatically when new cell types are added; full specifications are available in `models/threshold_config.json` in the source repository. A cross-domain contradiction checker automatically flags logical inconsistencies across the four domain assessments, including:

- high oncogenic potential with minimal network vulnerability
- high druggability with low clinical actionability
- an inhibition strategy conflicting with high tumor-suppressor likelihood
- a critical network node classified as easy to drug

The optional LLM layer (`USE_LLM_AGENTS=true`) adds natural-language rationales while preserving the same categorical output structure. An optional LLM narrative synthesis step (`USE_LLM_RECONCILIATION=true`) can synthesize cross-domain findings without introducing new assessments. All computational steps—network lookup, PageRank ranking, and pathway enrichment—are fully deterministic; LLM variability is isolated to the interpretation layer. The LLM abstraction layer supports Ollama (local or cloud), OpenAI, Anthropic, and any OpenAI-compatible API (e.g., Groq, Together, LM Studio), selected via the `LLM_PROVIDER` environment variable. The default Ollama path requires no API key, making the system self-contained.

The MCP server exposes thirteen tools, three browsable MCP Resources (cell-type catalogue, per-network summaries, and gene-lookup cards), and three MCP Prompt templates (single-gene deep dive, cancer biomarker panel, and cross-cell comparison) that scaffold common analysis workflows for new users, including gene validation (<100 ms), network queries (<50 ms) with optional ARACNe edge confidence filtering by mutual information score and bootstrap count, and reverse-direction master regulator analysis — given a set of differentially expressed genes, the `find_master_regulators` tool identifies which transcription factors most significantly drive that signature using Fisher's exact test enrichment against each regulon. New cell types can be incorporated through a documented pipeline covering data acquisition, ARACNe processing, cache generation, and automatic threshold computation; additional domain agents (e.g., epigenomics, metabolomics) could be integrated by extending the parallel analysis step. RegNetAgents runs on Windows, macOS, and Linux. LLM-generated domain interpretations should be treated as hypotheses for further investigation rather than validated conclusions. Interpretation quality depends on the completeness of the underlying ARACNe networks; genes absent from the network receive limited analysis. All deterministic outputs reported here are reproducible from the v1.0.1 source archived on Zenodo (DOI: 10.5281/zenodo.18500028), with dependencies pinned in the project configuration and an eleven-module test suite validating determinism, MCP integration, and agent-level behavior.

# Usage

RegNetAgents can be used programmatically via Python or through any MCP-compatible client (Claude Desktop, Cursor, Zed) via natural language (see installation guide):

```
Analyze the TP53 gene in epithelial cells
Compare MYC, TP53, and KRAS across different cell types
```

For reverse-direction analysis — identifying transcription factors that drive a gene signature — users supply a gene list in natural language:

```
These genes are upregulated in my RNA-seq experiment: TP53, CDKN1A, MDM2,
BAX, CCND1, MYC, EGFR, KRAS, CTNNB1, APC. Which transcription factors are
most likely driving this signature in epithelial cells?
```

## Example Application: Colorectal Cancer Gene Panel

The following output is from a five-gene colorectal cancer panel (MYC, CTNNB1, CCND1, TP53, KRAS in epithelial cells). The system identified 99 upstream regulators across all five genes and characterized each gene's network role, top upstream regulator by PageRank, and enriched pathway count:

```
Genes: MYC, CTNNB1, CCND1, TP53, KRAS | Cell Type: epithelial_cell

[OK] MYC:   hub_regulator | Regulators: 25 | Targets: 427 | Pathways: 58 | Top regulator: ID4    (PageRank: 0.622)
[OK] CTNNB1:hub_regulator | Regulators: 18 | Targets: 310 | Pathways:  2 | Top regulator: CHD2   (PageRank: 0.530)
[OK] CCND1: heavily_reg.  | Regulators: 42 | Targets:   0 | Pathways: 66 | Top regulator: ZBTB20 (PageRank: 0.600)
[OK] TP53:  hub_regulator | Regulators:  7 | Targets: 163 | Pathways: 16 | Top regulator: WWTR1  (PageRank: 0.473)
[OK] KRAS:  weakly_reg.   | Regulators:  7 | Targets:   0 | Pathways: 141| Top regulator: GPBP1  (PageRank: 0.609)
```

It can also be called programmatically without Claude Desktop. The repository includes `examples/quickstart.py`, a minimal single-gene entry point (2–5 seconds, no API keys required), and `demo_biomarker_analysis.py` for multi-gene panel analysis:

```python
import asyncio
from regnetagents_langgraph_workflow import RegNetAgentsWorkflow

async def main():
    workflow = RegNetAgentsWorkflow()
    result = await workflow.run_analysis(
        gene="TP53",
        cell_type="epithelial_cell",
        analysis_depth="comprehensive"
    )
    print(result["network_summary"])

asyncio.run(main())
```

The `analysis_depth` parameter accepts `"focused"` (network and pathways), `"basic"` (adds rule-based domain insights), or `"comprehensive"` (adds therapeutic prioritization and parallel domain analysis).

# Availability

RegNetAgents is available at <https://github.com/jab57/RegNetAgents> under the MIT License. The v1.0.1 release is archived on Zenodo (DOI: [10.5281/zenodo.18500028](https://doi.org/10.5281/zenodo.18500028)). Contribution guidelines and issue tracking are provided via GitHub.

# Acknowledgements

Network data were obtained from the GREmLN team's pre-computed ARACNe networks [@zhang2025gremln]. The author designed the system architecture, workflow logic, and reproducibility framework; Claude Code (Anthropic) assisted with implementation.

# References
