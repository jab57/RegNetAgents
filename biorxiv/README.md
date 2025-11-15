# bioRxiv Submission Package

> **Note:** This folder contains internal documentation for preparing the bioRxiv submission. These materials assume you have the RegNetAgents codebase for generating figures and testing. The software implementation is available from the corresponding author upon reasonable request for academic research purposes.

---

## Complete materials for submitting your preprint

---

## 📁 Files in This Folder

### **Core Submission Materials**

1. **`preprint_draft.md`** (in this folder)
   - Complete manuscript ready for bioRxiv (8,500 words)
   - Abstract, Introduction, Methods, Results, Discussion
   - 40 references included
   - Figure legends provided
   - **ACTION**: Edit author info (lines 8-10), then convert to PDF

2. **`generate_figures.py`** (located in: `../scripts/generate_figures.py`)
   - Python script to create figures from your data
   - Generates Figure 2 (biomarker panel) and Figure 3 (TP53 perturbation)
   - Also creates Tables 2 and 3 as CSV files
   - **ACTION**: Run with `python scripts/generate_figures.py`

### **Networking Materials (in this folder)**

3. **`one_page_summary.md`**
   - Concise 1-page research summary
   - Key results, validation, applications
   - Perfect for sharing with colleagues
   - **ACTION**: Convert to PDF to attach to emails

4. **`email_template.md`**
   - Ready-to-use email template for contacting 3 researchers
   - Instructions for finding researchers
   - Tips for customization
   - Response evaluation guide
   - **ACTION**: Customize and send to 3 people

### **Instructions (in this folder)**

5. **`submission_guide.md`**
   - Complete step-by-step bioRxiv submission instructions
   - Phase 1: Prepare materials (figures, PDF)
   - Phase 2: Submit to bioRxiv (15 minutes)
   - Phase 3: After posting (sharing strategy)
   - Troubleshooting common issues
   - **ACTION**: Follow this guide when ready to submit

6. **`README.md`** (this file)
   - Overview of all materials
   - Quick start guide

---

## 🚀 Quick Start (The Minimal Effort Plan)

### Week 1: Prepare and Submit (5 hours)

1. **Generate figures** (10 minutes)
   ```bash
   cd C:\Dev\RegNetAgents
   python scripts/generate_figures.py
   ```

2. **Review and edit draft** (30 minutes)
   - Open `preprint_draft.md` (in this folder)
   - Update lines 8-10 with your info
   - Skim through for any issues

3. **Create PDF** (15 minutes)
   - Use Pandoc, Word, or Google Docs
   - Embed generated figures
   - Save as `preprint.pdf`

4. **Submit to bioRxiv** (15-30 minutes)
   - Follow `submission_guide.md` in this folder
   - Create account at biorxiv.org
   - Upload PDF and fill out form

**Result**: Preprint posted with DOI in 1-3 days

---

### Week 2-3: Test Interest (1 hour)

5. **Prepare email materials** (15 minutes)
   - Convert `one_page_summary.md` to PDF
   - Add your bioRxiv URL to email template

6. **Send 3 emails** (30 minutes)
   - Follow `email_template.md` in this folder
   - Pick 3 researchers from suggestions
   - Send one email to each

7. **Wait 2 weeks** (0 hours)
   - No follow-up needed
   - Just wait to see if anyone responds

**Result**: Gauge interest in your work

---

### Week 4: Evaluate and Decide (30 minutes)

8. **Check responses**
   - 2-3 positive → Continue (optional Phase 4)
   - 0-1 responses → Stop here, move on

**Decision point**: You choose whether to invest more time

---

## 📊 What You'll Generate

When you run `generate_figures.py`, you'll create:

### Figures (for bioRxiv submission):
- `figure2_biomarker_panel.png` and `.pdf`
  - 4 panels: Regulatory architecture, therapeutic scores, biomarker types, roles
  - Based on your 5-gene CRC analysis

- `figure3_tp53_perturbation.png` and `.pdf`
  - 3 panels: Network diagram, PageRank ranking, degree centrality ranking
  - Shows 7 TP53 regulators with two alternative rankings (PageRank primary, degree secondary)

### Tables (for reference):
- `table2_biomarker_results.csv` and `.txt`
  - Summary of 5-gene analysis
  - Ready to copy into manuscript

- `table3_tp53_perturbation.csv` and `.txt`
  - TP53 perturbation results
  - 7 regulators with impact scores

---

## ✅ Pre-Submission Checklist

Before you submit to bioRxiv:

- [ ] Run `python scripts/generate_figures.py` successfully
- [ ] Edit `preprint_draft.md` with your name/email/institution
- [ ] Review abstract and make sure you're comfortable with all claims
- [ ] Create PDF with embedded figures
- [ ] Check PDF looks good (no formatting issues)
- [ ] Create bioRxiv account
- [ ] Have abstract, keywords, and data availability statement ready (all in draft)

---

## 📧 Email Checklist

After preprint is posted:

- [ ] Convert `one_page_summary.md` to PDF
- [ ] Add your bioRxiv URL to email template
- [ ] Identify 3 researchers to contact
- [ ] Customize email template for each person
- [ ] Send 3 emails
- [ ] Set calendar reminder for 2 weeks to check responses

---

## 🎯 Success Criteria

You'll know this was worth it if:

✅ **Minimum success** (guaranteed):
- You have a citable preprint with DOI
- Work is documented and discoverable
- Priority is established for your innovation

✅ **Ideal success** (if Phase 3 goes well):
- 2-3 researchers respond positively
- Someone asks about collaboration
- Validation that work is interesting to others

---

## ⏱️ Time Commitment Summary

| Phase | Time | When |
|-------|------|------|
| Prepare materials | 2 hours | Week 1 |
| Submit to bioRxiv | 30 min | Week 1 |
| Send 3 emails | 30 min | Week 2 |
| Evaluate responses | 30 min | Week 4 |
| **TOTAL** | **~4 hours** | **1 month** |

**Optional Phase 4** (if interest is high): 2-10 more hours

---

## 🛠️ Requirements

### Software Needed:

**To generate figures**:
- Python 3.8+
- matplotlib
- pandas
- numpy

Install with:
```bash
pip install matplotlib pandas numpy
```

**To create PDF** (pick one):
- Pandoc (free: https://pandoc.org)
- Microsoft Word or Google Docs
- Overleaf (free LaTeX: https://overleaf.com)

---

## 📁 Project Structure

```
RegNetAgents/
├── README.md
├── requirements.txt
├── regnetagents_langgraph_mcp_server.py (main server)
├── regnetagents_langgraph_workflow.py (main workflow)
├── gene_id_mapper.py
├── complete_gene_service.py
│
├── scripts/ (utility scripts)
│   └── generate_figures.py
│
├── results/ (analysis outputs)
│   ├── biomarker_results.json
│   └── tp53_perturbation_results.json
│
├── docs/ (documentation)
├── tests/ (test files)
├── models/ (network data)
├── agents/ (agent definitions)
│
└── biorxiv/ (THIS FOLDER)
    ├── README.md (this file)
    ├── preprint_draft.md (manuscript - 8,500 words)
    ├── submission_guide.md (step-by-step instructions)
    ├── email_template.md (networking guide)
    └── one_page_summary.md (share with colleagues)
```

---

## 🆘 Getting Help

### During figure generation:
- Check Python is installed: `python --version`
- Check packages: `pip list`
- Error messages: Usually about missing packages → `pip install <package>`

### During PDF creation:
- Pandoc issues: Check installation with `pandoc --version`
- Word/Docs: Just copy/paste markdown → export as PDF
- Can't embed figures: Use a PDF editor or insert manually

### During submission:
- Read `submission_guide.md` carefully
- bioRxiv FAQ: https://www.biorxiv.org/about/FAQ
- Email bioRxiv: submit@biorxiv.org

---

## 🎉 You're Ready!

Everything you need is organized:
- ✅ Complete manuscript (this folder)
- ✅ Figure generation code (parent folder)
- ✅ Submission instructions (this folder)
- ✅ Networking templates (this folder)
- ✅ Clear timeline

**Next step**: Open `submission_guide.md` and start Phase 1 when you're ready.

**Remember**: This is low-risk, low-commitment. You can stop at any phase if it's not working for you.

Good luck! 🚀
