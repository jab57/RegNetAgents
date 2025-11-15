# RegNetAgents Documentation

## Complete documentation for the RegNetAgents framework

---

## 📚 Documentation Files

### **Setup & Getting Started**

1. **REGNETAGENTS_MCP_SETUP.md**
   - How to install and configure the MCP server
   - Claude Desktop integration
   - Quick start guide

### **Architecture & Design**

2. **REGNETAGENTS_Analysis_Pipeline.md**
   - Complete workflow architecture
   - Multi-agent system design
   - LangGraph orchestration details

3. **GENE_MAPPING_ARCHITECTURE.md**
   - Gene ID mapping system
   - Symbol ↔ Ensembl conversion
   - Caching strategy

### **Data Sources & Processing**

4. **DATA_SOURCES.md** ⭐ NEW
   - Complete documentation of the 10 current cell types
   - Source datasets from CellxGene Portal
   - Cell type characteristics and research applications
   - How to find and download similar datasets
   - Data quality requirements

5. **END_TO_END_DATA_PIPELINE.md** ⭐ NEW
   - Complete pipeline from CellxGene to ready-to-use caches
   - Step-by-step walkthrough of all processing stages
   - Quality control and preprocessing details
   - ARACNe network inference process
   - Cache building and PageRank pre-computation
   - Performance benchmarks and troubleshooting

### **Adding Cell Types**

6. **ADDING_NEW_CELL_TYPES.md** ⭐ CONSOLIDATED
   - Complete guide to add new cell types
   - Data acquisition from CellxGene Portal
   - Recommended datasets for 5 future cell types (hepatocytes, cardiomyocytes, neurons, fibroblasts, endothelial cells)
   - Quality control criteria and requirements
   - HPC setup and ARACNe pipeline execution
   - Validation, testing, and integration
   - Timeline and resource planning

### **Integration & Features**

7. **REACTOME_INTEGRATION_SUMMARY.md**
   - Reactome API integration
   - Pathway enrichment analysis
   - Statistical validation methods

### **Conference & Dissemination**

8. **CONFERENCE_SUBMISSION_GUIDE_2026.md**
   - List of relevant conferences
   - Submission strategies
   - Abstract writing tips
   - Top recommendations: ISMB, Single Cell Genomics, ICSB

9. **REGNETAGENTS_CONFERENCE_POSTER.md**
   - Poster content and design
   - Key results and figures
   - Technical details for presentation

10. **REGNETAGENTS_CONFERENCE_POSTER.pptx**
    - PowerPoint version of poster
    - Ready for printing/presenting

---

## 📖 Quick Reference Guide

### For Users:
- **Getting Started**: Start with `REGNETAGENTS_MCP_SETUP.md`
- **Understanding the System**: Read `REGNETAGENTS_Analysis_Pipeline.md`

### For Developers:
- **Architecture**: `REGNETAGENTS_Analysis_Pipeline.md` + `GENE_MAPPING_ARCHITECTURE.md`
- **Adding Features**: See workflow code + architecture docs

### For Researchers:
- **Publishing**: Check `CONFERENCE_SUBMISSION_GUIDE_2026.md`
- **Poster**: Use `REGNETAGENTS_CONFERENCE_POSTER.md` or `.pptx`

### For Data Scientists:
- **Adding Cell Types**: `ADDING_NEW_CELL_TYPES.md` (comprehensive guide)
- **Data Pipeline**: `END_TO_END_DATA_PIPELINE.md`

---

## 📊 Documentation by Topic

### **Installation & Setup**
- REGNETAGENTS_MCP_SETUP.md

### **System Architecture**
- REGNETAGENTS_Analysis_Pipeline.md
- GENE_MAPPING_ARCHITECTURE.md

### **Data Sources & Pipeline** ⭐
- DATA_SOURCES.md (comprehensive source documentation)
- END_TO_END_DATA_PIPELINE.md (complete processing pipeline)

### **Extending the System**
- ADDING_NEW_CELL_TYPES.md (complete guide for adding new cell types)
- REACTOME_INTEGRATION_SUMMARY.md (pathway enrichment)

### **Publication & Dissemination**
- CONFERENCE_SUBMISSION_GUIDE_2026.md
- REGNETAGENTS_CONFERENCE_POSTER.md
- REGNETAGENTS_CONFERENCE_POSTER.pptx

---

## 🔗 Related Folders

- **`../biorxiv/`** - bioRxiv preprint materials
- **`../scripts/`** - Utility scripts referenced in docs
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
4. **CONFERENCE_SUBMISSION_GUIDE_2026.md** - Share your work

For adding new cell types:
1. **DATA_SOURCES.md** - Learn about current data
2. **END_TO_END_DATA_PIPELINE.md** - Complete processing walkthrough
3. **ADDING_NEW_CELL_TYPES.md** - Comprehensive activation guide

---

## 📧 Need Help?

- Check the specific documentation file for your topic
- Review the main `../README.md` for project overview
- See `../biorxiv/README.md` for publication guidance

---

Last Updated: 2025-11-08
