# RegNetAgents Documentation

## Complete documentation for the RegNetAgents framework

---

## 📚 Documentation Files

### **Setup & Getting Started**

1. **REGNETAGENTS_MCP_SETUP.md**
   - How to install and configure the MCP server
   - Claude Desktop integration
   - Quick start guide

2. **CLAUDE_DESKTOP_USAGE.md**
   - Detailed Claude Desktop usage guide
   - Example queries and workflows

### **Architecture & Design**

3. **REGNETAGENTS_Analysis_Pipeline.md**
   - Complete workflow architecture
   - Multi-agent system design
   - LangGraph orchestration details

4. **GENE_MAPPING_ARCHITECTURE.md**
   - Gene ID mapping system
   - Symbol ↔ Ensembl conversion
   - Caching strategy

### **Data Sources & Processing**

5. **DATA_SOURCES.md**
   - Complete documentation of the 10 current cell types
   - Source datasets from CellxGene Portal
   - Cell type characteristics and research applications
   - How to find and download similar datasets
   - Data quality requirements

6. **END_TO_END_DATA_PIPELINE.md**
   - Complete pipeline from CellxGene to ready-to-use caches
   - Step-by-step walkthrough of all processing stages
   - Quality control and preprocessing details
   - ARACNe network inference process
   - Cache building and PageRank pre-computation
   - Performance benchmarks and troubleshooting

### **Adding Cell Types**

7. **ADDING_NEW_CELL_TYPES.md**
   - Complete guide to add new cell types
   - Data acquisition from CellxGene Portal
   - Quality control criteria and requirements
   - HPC setup and ARACNe pipeline execution
   - Validation, testing, and integration

### **Integration & Features**

8. **REACTOME_INTEGRATION_SUMMARY.md**
   - Reactome API integration
   - Pathway enrichment analysis
   - Statistical validation methods

---

## 📖 Quick Reference Guide

### For Users:
- **Getting Started**: Start with `REGNETAGENTS_MCP_SETUP.md`
- **Understanding the System**: Read `REGNETAGENTS_Analysis_Pipeline.md`

### For Developers:
- **Architecture**: `REGNETAGENTS_Analysis_Pipeline.md` + `GENE_MAPPING_ARCHITECTURE.md`
- **Adding Features**: See workflow code + architecture docs

### For Data Scientists:
- **Adding Cell Types**: `ADDING_NEW_CELL_TYPES.md` (comprehensive guide)
- **Data Pipeline**: `END_TO_END_DATA_PIPELINE.md`

---

## 📊 Documentation by Topic

### **Installation & Setup**
- REGNETAGENTS_MCP_SETUP.md
- CLAUDE_DESKTOP_USAGE.md

### **System Architecture**
- REGNETAGENTS_Analysis_Pipeline.md
- GENE_MAPPING_ARCHITECTURE.md

### **Data Sources & Pipeline**
- DATA_SOURCES.md (comprehensive source documentation)
- END_TO_END_DATA_PIPELINE.md (complete processing pipeline)

### **Extending the System**
- ADDING_NEW_CELL_TYPES.md (complete guide for adding new cell types)
- REACTOME_INTEGRATION_SUMMARY.md (pathway enrichment)

---

## 🔗 Related Folders

- **`../scripts/`** - Utility scripts (build_network_cache.py, build_gene_annotation_database.py)
- **`../tests/`** - Test files for validation

---

## 📝 Documentation Standards

All documentation follows these conventions:
- **Markdown format** for readability
- **Code blocks** for commands and examples
- **Step-by-step instructions** where applicable
- **Clear section headers** for navigation

---

## 🎯 Most Important Documents

If you're short on time, read these first:

1. **REGNETAGENTS_MCP_SETUP.md** - Get the system running
2. **DATA_SOURCES.md** - Understand the data (10 current cell types)
3. **REGNETAGENTS_Analysis_Pipeline.md** - How the system works

For adding new cell types:
1. **DATA_SOURCES.md** - Learn about current data
2. **END_TO_END_DATA_PIPELINE.md** - Complete processing walkthrough
3. **ADDING_NEW_CELL_TYPES.md** - Comprehensive activation guide

---

## 📧 Need Help?

- Check the specific documentation file for your topic
- Review the main `../README.md` for project overview

---

Last Updated: 2026-01-31
