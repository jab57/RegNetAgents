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
    orcid: 0000-0000-0000-0000
    affiliation: 1
affiliations:
  - name: Bird AI Solutions
    index: 1
date: 2 February 2025
bibliography: paper.bib
---

# Summary

RegNetAgents is a Python-based multi-agent system that automates gene regulatory network analysis through a hybrid architecture combining large language model (LLM) reasoning with deterministic rule-based computation. Built on the LangGraph orchestration framework [@langgraph], the system coordinates four specialized domain agents (cancer biology, drug discovery, clinical relevance, systems biology) that analyze gene regulatory relationships derived from ARACNe-inferred networks [@margolin2006aracne; @lachmann2016aracne]. The framework integrates with Claude Desktop via the Model Context Protocol (MCP) [@mcp], enabling researchers to perform sophisticated network analyses through natural language queries.

The system processes pre-computed ARACNe networks from the GREmLN project [@zhang2025gremln], covering 10 cell types derived from over 500,000 single cells in the CELLxGENE Data Portal [@megill2021cellxgene]. For each query gene, RegNetAgents identifies upstream regulators and downstream targets, performs Reactome pathway enrichment [@gillespie2022reactome], calculates network centrality metrics (PageRank, degree centrality) for therapeutic target prioritization [@mora2021effects], and synthesizes multi-domain biological interpretations.

Key features include:

- **Hybrid LLM+deterministic architecture**: Graceful fallback from LLM-powered insights to rule-based logic ensures reliable execution regardless of LLM availability
- **Parallel multi-agent execution**: Four domain agents operate concurrently, reducing analysis time from hours to seconds
- **Therapeutic target prioritization**: PageRank-based ranking of upstream regulators identifies candidates for experimental validation
- **Conversational interface**: MCP integration enables natural language interaction through Claude Desktop

# Statement of Need

Gene regulatory network analysis is fundamental to understanding cellular function, disease mechanisms, and therapeutic opportunities [@califano2017recurrent; @sonawane2017understanding]. However, traditional workflows require researchers to manually query multiple databases (STRING, BioGRID, Reactome, Enrichr), export and integrate data across platforms, and synthesize findings across specialized domains. This fragmented process typically requires 2-4 hours per gene and demands both programming expertise and deep domain knowledge.

RegNetAgents addresses these barriers by automating the complete analytical workflow through a conversational interface. The system transforms multi-hour, multi-tool pipelines into second-scale execution: a five-gene colorectal cancer biomarker panel (MYC, CTNNB1, CCND1, TP53, KRAS) is analyzed in 15-62 seconds compared to an estimated 12+ hours manually. This democratizes access to sophisticated network analysis methods for researchers without computational backgrounds.

The framework is designed for hypothesis generation and experimental prioritization rather than prediction of molecular outcomes. Validation against a well-characterized cancer biomarker panel demonstrates the system recapitulates literature-confirmed regulatory relationships (e.g., WWTR1/YAP1 regulation of TP53 [@strano2005yap; @levy2008yap1]) while generating testable hypotheses for experimental follow-up.

No existing tool provides this combination of multi-agent LLM orchestration, deterministic fallback guarantees, and conversational access to gene regulatory network analysis. Related systems like AutoGPT [@autogpt], CAMEL [@li2023camel], and MetaGPT [@hong2023metagpt] provide general-purpose agent frameworks but lack domain-specific scientific workflow integration and reliability guarantees for reproducible research.

# Implementation

RegNetAgents is implemented in Python 3.8+ using LangGraph for workflow orchestration, NetworkX [@hagberg2008networkx] for graph algorithms, and the MCP SDK for Claude Desktop integration. The directed acyclic graph (DAG) workflow comprises:

1. **Initialization**: Gene symbol validation and Ensembl ID resolution
2. **Network modeling**: Regulator/target identification from ARACNe networks with role classification
3. **Therapeutic prioritization**: PageRank-based ranking of upstream regulators (triggered for genes with >5 regulators)
4. **Pathway enrichment**: Reactome API queries with FDR correction
5. **Domain analysis**: Parallel execution of four specialized agents
6. **Integration**: Synthesis of findings into structured JSON reports

Domain agents support dual-mode operation: LLM-powered (Ollama with llama3.1:8b) for nuanced biological interpretation, or rule-based for deterministic execution. All network centrality metrics are computed algorithmically regardless of mode, ensuring reproducible quantitative outputs.

# Acknowledgements

This project was developed with assistance from Claude Code (Anthropic). Network data were obtained from the GREmLN team's pre-computed ARACNe networks [@zhang2025gremln].

# References
