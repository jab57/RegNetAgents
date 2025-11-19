# RegNetAgents: Multi-Agent AI for Gene Regulatory Network Analysis

**One-Page Research Summary**

---

## The Problem

Gene regulatory network analysis currently requires:
- ⏱️ **Multiple hours** of manual work per gene
- 🔀 **Multiple tools**: STRING, BioGRID, Enrichr, DAVID, manual literature curation
- 📊 **Sequential processing**: Analyzing 5 genes = 5× the time
- 🧠 **Expert knowledge**: Separate analysis for cancer, drug, clinical, systems biology

**Result**: Slow, fragmented, inaccessible to non-computational researchers

---

## Our Solution: RegNetAgents

**LLM-powered multi-agent AI framework** that automates comprehensive gene regulatory analysis:

### Key Innovations

1. **⚡ Fast Analysis**
   - Rule-based mode: **15.49 seconds** for 5 genes (99 regulators analyzed)
   - LLM-powered mode: **~62 seconds** for 5 genes with scientific rationales
   - 115-480× faster than manual workflows (hours → seconds)

2. **🤖 LLM-Powered Multi-Agent System**
   - 4 specialized domain agents (cancer, drug, clinical, systems biology)
   - Local LLM inference (Ollama/llama3.1:8b) generates scientific insights with rationales
   - Execute in parallel for integrated analysis
   - Graceful fallback to rule-based heuristics for reliability

3. **🎯 Perturbation Analysis**
   - Simulates inhibiting upstream regulators (topology-based, not expression prediction)
   - Ranks therapeutic targets using PageRank and out-degree centrality for experimental validation
   - Identifies novel drug targets through network analysis

4. **💬 Conversational Interface**
   - Natural language queries via Claude Desktop
   - No programming required
   - Accessible to experimental biologists

---

## Demonstration: Colorectal Cancer Biomarkers

**Case Study**: 5-gene panel (MYC, CTNNB1, CCND1, TP53, KRAS)

**Results**:
- ✅ **Framework validation** - **100% concordance** with published literature across five genes
- ✅ **Automated classification** - diagnostic, prognostic, predictive biomarkers
- ✅ **Perturbation analysis** - all 99 regulators analyzed across 5 genes
- ✅ **Validated predictions** - successfully identified experimentally validated TP53 regulators (WWTR1, YAP1, CHD4 from Hippo pathway)
- ✅ **15.49 seconds** - complete network and perturbation analysis (all 99 regulators, rule-based)
- ✅ **~62 seconds** - comprehensive analysis with LLM-powered domain insights

### TP53 Perturbation Analysis Example

**Prioritized 7 candidate regulators** for experimental validation:
- **WWTR1**: Top-ranked by PageRank (0.473) - Literature confirms TP53 interaction ✓ (Hippo pathway)
- **CHD4**: Literature-confirmed TP53 regulator ✓ (regulates TP53 acetylation)
- **YAP1**: Literature-confirmed regulator ✓ (Hippo pathway effector)
- **RBPMS**: Testable hypothesis (highest out-degree centrality: 403 downstream targets)

All candidates ranked by PageRank (primary) and out-degree centrality (secondary) for experimental validation prioritization.

**Note**: Perturbation analysis uses network topology to prioritize regulators for experimental validation. It does not predict gene expression levels or dynamic responses. Framework designed for hypothesis generation, not therapeutic claims.

---

## Technology

**Network Data**: Pre-computed regulatory networks from GREmLN foundation model (11M scRNA-seq profiles, 162 cell types)
**Network Inference**: ARACNe algorithm via GREmLN team preprocessing (statistical mutual information method)
**Architecture**: LangGraph workflow orchestration + Model Context Protocol
**LLM Engine**: Ollama (local inference, llama3.1:8b) with graceful fallback
**Data Sources**: GREmLN/CELLxGENE Data Portal (10 cell types), Reactome Pathway Database
**Processing**: GREmLN team (Zhang et al. 2025, CZ Biohub NY / Columbia University)

---

## Key Features

| Feature | Traditional | RegNetAgents |
|---------|------------|---------------|
| **Speed** | Hours | Seconds (15-62 sec for 5 genes) |
| **Integration** | Manual | Automated |
| **Domains** | Siloed | 4 LLM-powered parallel agents |
| **Insights** | Manual literature | AI-generated with rationales |
| **Interface** | Web forms | Conversational |
| **Perturbation** | Manual | Automated |
| **Cell types** | Repeat workflow | Instant (<0.01 sec) |

---

## Applications

✅ **Hypothesis Generation**: Prioritize candidate regulators for experimental validation
✅ **Experimental Design**: Identify testable hypotheses from network topology
✅ **Disease Mechanisms**: Cross-cell-type regulatory analysis
✅ **Workflow Automation**: Replace multi-hour manual analysis with second-scale queries

---

## Availability

**Preprint**: bioRxiv [DOI to be added]
**Data**:
- GREmLN preprocessed networks: [Quickstart Tutorial](https://virtualcellmodels.cziscience.com/quickstart/gremln-quickstart)
- Underlying data: CELLxGENE Census 2024-07-01 + Reactome pathways

**Development Note**: This framework was developed with significant assistance from Claude Code (Anthropic's AI coding agent), demonstrating how AI-assisted development tools can enable researchers from diverse backgrounds to tackle complex computational biology challenges.

---

## Contact

Jose A. Bird, PhD
Bird AI Solutions
jbird@birdaisolutions.com

---

## Key Results at a Glance

**Performance**:
- Rule-based mode (5 genes, 99 regulators): 15.49 seconds
- LLM-powered mode (5 genes with rationales): ~62 seconds
- Single gene (LLM mode): ~15 seconds
- Single gene (rule-based): 0.60 seconds
- Cross-cell (10 types): <0.01 seconds
- 115-480× faster than manual literature review

**Framework Validation**:
- Regulatory patterns align with published literature
- Top-ranked TP53 regulators include experimentally validated Hippo pathway effectors
- Network topology-based ranking recapitulates known biology and generates testable hypotheses

**Innovation**:
- Multi-agent framework with LLM-powered domain analysis
- Local inference (Ollama) with graceful fallback to rules
- AI-generated scientific rationales for all domain insights
- Automated perturbation analysis for candidate regulator prioritization
- Conversational interface making network biology accessible via Model Context Protocol
- Hypothesis generation tool for experimental prioritization

---

**Transform hours of manual analysis into seconds of automated insights.**
