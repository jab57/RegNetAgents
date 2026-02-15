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

RegNetAgents is a Python package that automates gene regulatory network analysis using a multi-agent architecture built on LangGraph [@langgraph]. The software coordinates four specialized domain agents (cancer biology, drug discovery, clinical relevance, systems biology) to analyze regulatory relationships from ARACNe-inferred networks [@margolin2006aracne; @lachmann2016aracne]. It integrates with Claude Desktop via the Model Context Protocol (MCP) [@mcp], enabling natural language queries for network analysis. RegNetAgents is intended for computational biologists, bioinformaticians, and researchers who need automated GRN analysis without writing code.

The architecture separates workflow logic from LLM interpretation: a fixed directed acyclic graph (DAG) controls execution order and data flow, while LLMs provide domain-specific synthesis of results. This separation ensures that identical inputs produce identical computational outputs, with LLM variability isolated to the interpretation layer.

Key features:

- **Hybrid architecture**: Combines LLM-powered interpretation with deterministic rule-based fallback for reliable execution
- **Parallel agents**: Four domain agents run concurrently for fast analysis
- **Therapeutic prioritization**: Deterministic PageRank and centrality-based ranking of upstream regulators [@mora2021effects]
- **Pathway enrichment**: Automated Reactome API queries with FDR correction [@gillespie2022reactome]

# Statement of Need

Researchers analyzing gene regulatory networks typically query multiple databases manually (STRING, BioGRID, Reactome), export data across platforms, and synthesize findings—a process requiring hours per gene and programming expertise. Established tools address individual components of this workflow: Cytoscape [@shannon2003cytoscape] provides interactive network visualization, pySCENIC [@aibar2017scenic] and CellOracle [@kamimoto2023dissecting] infer regulatory networks from expression data, and STRING [@szklarczyk2023string] catalogs known protein interactions. However, these tools operate on separate steps that must be manually connected, require programming expertise, and produce outputs that need expert interpretation across multiple biological contexts. These tools focus on network inference or visualization, whereas RegNetAgents operates downstream, automating multi-perspective interpretation of pre-computed networks. **No existing tool in this landscape provides an integrated, conversational workflow that combines GRN analysis with deterministic reproducibility.**

RegNetAgents addresses this gap by combining a natural-language interface with a reproducible, multi-agent workflow that performs regulator/target identification, therapeutic prioritization, pathway enrichment, and domain-specific interpretation in a single automated pipeline. Unlike network inference tools (pySCENIC, CellOracle) that construct networks from raw data, RegNetAgents operates downstream—analyzing pre-computed networks to answer biological questions without requiring users to write code. In a five-gene colorectal cancer case study, the system analyzed 99 upstream regulators in 15–62 seconds. The equivalent manual workflow—querying regulators and targets across 10 cell types, performing pathway enrichment for each, and synthesizing findings across four biological perspectives per gene—requires an estimated 8–16 hours, a speedup of three orders of magnitude.

The software processes pre-computed ARACNe networks from the GREmLN project [@zhang2025gremln], covering 10 cell types from the CELLxGENE Data Portal [@megill2021cellxgene]. Unlike general-purpose agent frameworks [@autogpt; @li2023camel; @hong2023metagpt], RegNetAgents provides domain-specific scientific workflow integration with deterministic fallback guarantees for reproducible research.

# Implementation

![RegNetAgents multi-agent architecture.](figure1_architecture.png)

RegNetAgents requires Python 3.10+ and uses NetworkX [@hagberg2008networkx] for graph algorithms. The workflow is implemented as a LangGraph DAG, which enforces explicit stepwise execution and ensures reproducibility across runs. This design intentionally reduces agent autonomy—agents cannot dynamically alter the workflow—in exchange for deterministic, auditable execution. Three design trade-offs follow from this choice: the DAG constraints limit agent autonomy but guarantee reproducibility and auditability; the dual-mode architecture (LLM + rule-based) trades richer interpretation for guaranteed determinism when needed; and restricting input to pre-computed ARACNe networks prioritizes validated, curated data over generality. The pipeline includes gene validation, network lookup, therapeutic prioritization (for genes with >5 regulators), pathway enrichment, parallel domain analysis, and report generation.

Therapeutic prioritization uses NetworkX's deterministic PageRank and degree-based centrality metrics to rank upstream regulators by influence within the ARACNe-derived subnetwork. Domain agents support dual-mode operation—LLM-powered or rule-based—with the rule-based fallback ensuring reproducible outputs even when LLM variability is present. All computational steps—network lookup, PageRank ranking, and pathway enrichment—are fully deterministic; only the natural-language interpretation layer varies between runs. The LLM layer uses Ollama, which supports local or cloud deployment and runs a range of open-weight models (e.g., Llama, Mistral); the modular architecture allows alternative LLM backends to be substituted.

The MCP server exposes eleven tools including a lightweight gene validation tool that checks names against the network in under 100 ms and returns fuzzy-matched suggestions for misspelled symbols, and a network query tool that answers structural questions (top regulators, gene neighbors, network statistics) in under 50 ms. The architecture is modular, allowing users to add new domain agents or replace components without modifying the core DAG. The package includes installation verification, a test suite with eight modules validating MCP integration, workflow determinism, and agent-level behavior, and documentation covering setup, API usage, and tutorials. A minimal usage example in the README enables users to run a complete GRN analysis with a single command. RegNetAgents runs on Windows, macOS, and Linux. Current network coverage spans 10 cell types from the GREmLN project; additional cell types can be incorporated through a documented pipeline covering data acquisition, ARACNe processing, and cache generation. LLM-generated domain interpretations should be treated as hypotheses for further investigation rather than validated conclusions.

# Research Application

To illustrate practical usage, a five-gene colorectal cancer case study demonstrates that RegNetAgents completes a full analysis pipeline—from gene validation through multi-perspective domain synthesis—in 15–62 seconds per gene. For a representative query on TP53, the system identifies 7 upstream regulators and 163 downstream targets across 10 cell types, performs PageRank-based therapeutic prioritization, runs pathway enrichment against Reactome, and synthesizes findings from four domain perspectives in a single automated run.

In rule-based mode, all outputs—regulator lists, rankings, enrichment results, and domain reports—are identical across runs, supporting reproducible computational workflows. The modular architecture also supports extensibility: new cell types can be added through the documented data pipeline without modifying the core analysis logic.

# Availability

RegNetAgents is available at <https://github.com/jab57/RegNetAgents> under the MIT License. Contribution guidelines and issue tracking are provided via GitHub.

# Acknowledgements

Network data were obtained from the GREmLN team's pre-computed ARACNe networks [@zhang2025gremln]. The author designed the system architecture, workflow logic, and reproducibility framework; Claude Code (Anthropic) assisted with implementation.

# References
