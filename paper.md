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
date: 2 February 2025
bibliography: paper.bib
---

# Summary

RegNetAgents is a Python package that automates gene regulatory network analysis using a multi-agent architecture built on LangGraph [@langgraph]. The software coordinates four specialized domain agents (cancer biology, drug discovery, clinical relevance, systems biology) to analyze regulatory relationships from ARACNe-inferred networks [@margolin2006aracne; @lachmann2016aracne]. It integrates with Claude Desktop via the Model Context Protocol (MCP) [@mcp], enabling natural language queries for network analysis. RegNetAgents is intended for computational biologists, bioinformaticians, and researchers who need automated GRN analysis without writing code.

Key features:

- **Hybrid architecture**: Combines LLM-powered interpretation with deterministic rule-based fallback for reliable execution
- **Parallel agents**: Four domain agents run concurrently for fast analysis
- **Therapeutic prioritization**: PageRank-based ranking of upstream regulators [@mora2021effects]
- **Pathway enrichment**: Automated Reactome API queries with FDR correction [@gillespie2022reactome]

# Statement of Need

Researchers analyzing gene regulatory networks typically query multiple databases manually (STRING, BioGRID, Reactome), export data across platforms, and synthesize findings—a process requiring hours per gene and programming expertise. No existing tool provides an integrated, conversational workflow for GRN analysis with deterministic reproducibility. RegNetAgents addresses this gap by automating the complete workflow through a conversational interface, reducing analysis time from hours to seconds.

The software processes pre-computed ARACNe networks from the GREmLN project [@zhang2025gremln], covering 10 cell types from the CELLxGENE Data Portal [@megill2021cellxgene]. Unlike general-purpose agent frameworks [@autogpt; @li2023camel; @hong2023metagpt], RegNetAgents provides domain-specific scientific workflow integration with deterministic fallback guarantees for reproducible research.

# Implementation

![RegNetAgents multi-agent architecture showing the LangGraph-orchestrated workflow from user query through domain analysis to structured output.](figure1_architecture.png)

RegNetAgents requires Python 3.10+ and uses NetworkX [@hagberg2008networkx] for graph algorithms. The workflow is implemented as a LangGraph DAG, ensuring reproducible, stepwise execution: gene validation, network lookup, therapeutic prioritization (for genes with >5 regulators), pathway enrichment, parallel domain analysis, and report generation. Domain agents support dual-mode operation—LLM-powered (Ollama) or rule-based—with all centrality metrics computed deterministically. The package includes installation verification, a test suite with 10 test modules, and documentation covering setup, API usage, and tutorials. RegNetAgents runs on Windows, macOS, and Linux.

# Availability

RegNetAgents is available at <https://github.com/jab57/RegNetAgents> under the MIT License.

# Acknowledgements

Network data were obtained from the GREmLN team's pre-computed ARACNe networks [@zhang2025gremln]. This project was developed with assistance from Claude Code (Anthropic).

# References
