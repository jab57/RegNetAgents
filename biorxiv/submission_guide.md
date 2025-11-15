# bioRxiv Submission Guide

## Complete Step-by-Step Instructions

This guide will walk you through submitting your preprint to bioRxiv in ~20 minutes.

---

## 🚀 Quick Start - What to Do RIGHT NOW

Most prep work is already complete! You just need to:

1. **Convert to PDF** (2 min): Open `regnetagents_preprint.docx` → Save As PDF
2. **Review PDF** (3 min): Make sure everything looks good
3. **Go to Phase 2** (15 min): Follow bioRxiv submission steps below

**That's it!** All figures, tables, and manuscript content are ready.

---

## Phase 1: Prepare Your Materials (Before Submission)

### ✅ ALREADY COMPLETE

The following materials are ready in `C:\Dev\RegNetAgents\biorxiv\`:

**Figures (3 complete sets)**:
- ✅ `figure1_architecture.png` and `.pdf`
- ✅ `figure2_biomarker_panel.png` and `.pdf`
- ✅ `figure3_tp53_perturbation.png` and `.pdf`

**Tables (2 complete sets)**:
- ✅ `table2_biomarker_results.csv` and `.txt`
- ✅ `table3_tp53_perturbation.csv` and `.txt`

**Manuscript**:
- ✅ `preprint_draft.md` - Updated with author information
- ✅ `regnetagents_preprint.docx` - Word document ready for conversion

---

### Step 1: Create PDF (~2 minutes)

You already have `regnetagents_preprint.docx` ready to convert.

**Method: Convert Word to PDF**

1. Open `C:\Dev\RegNetAgents\biorxiv\regnetagents_preprint.docx` in Microsoft Word
2. Click **File** → **Save As**
3. Choose location: Same folder (`biorxiv`)
4. File name: `regnetagents_preprint.pdf`
5. Save as type: **PDF**
6. Click **Save**

**Alternative: If you don't have Word installed**

Use Word Online (free):
1. Go to https://www.office.com
2. Sign in with Microsoft account (or create free account)
3. Upload `regnetagents_preprint.docx`
4. Click **File** → **Download as PDF**

**Final check** - Open the PDF and verify:
- ✅ All text is readable
- ✅ Title and headers are properly formatted
- ✅ References are included
- ✅ Your name shows: "Jose A. Bird, PhD"
- ✅ Your affiliation shows: "Founder & Data Science Consultant, Bird AI Solutions"

**Note**: The Word document doesn't include embedded figures. You'll upload figures separately to bioRxiv as supplementary files (instructions in Phase 2).

---

## Phase 2: Submit to bioRxiv (~15 minutes)

### Step 1: Create Account

1. Go to https://www.biorxiv.org
2. Click "Submit a Manuscript" (top right)
3. Create account with your email
4. Verify email address

---

### Step 2: Start Submission

1. Log in to bioRxiv
2. Click "Submit a Manuscript"
3. Click "New Submission"

---

### Step 3: Fill Out Submission Form

#### Part 1: Manuscript Type
- Select: **"New Results"**

#### Part 2: Subject Category
- Primary: **"Bioinformatics"**
- Secondary (optional): **"Systems Biology"**, **"Cancer Biology"**

#### Part 3: Upload Manuscript
- Click "Choose File"
- Select `regnetagents_preprint.pdf` from your biorxiv folder
- Wait for upload (may take 1-2 minutes)

#### Part 3b: Upload Figures (Supplementary Files)
- Look for "Supplementary Files" or "Additional Files" section
- Upload the following figures:
  - `figure1_architecture.pdf`
  - `figure2_biomarker_panel.pdf`
  - `figure3_tp53_perturbation.pdf`
- Label them as "Figure 1", "Figure 2", "Figure 3" respectively

#### Part 4: Manuscript Information

**Title** (copy from your PDF):
```
RegNetAgents: A Multi-Agent AI Framework for Automated Gene Regulatory Network Analysis
```

**Running Title**:
```
Multi-Agent Gene Regulatory Network Analysis
```

**Abstract** (copy from the updated preprint draft):
```
Gene regulatory network analysis is essential for understanding disease mechanisms, identifying biomarkers, and prioritizing therapeutic targets. However, traditional workflows require hours of manual effort across multiple databases and tools, with sequential processing limiting scalability. We present RegNetAgents, an LLM-powered multi-agent artificial intelligence framework that automates comprehensive gene regulatory analysis through intelligent workflow orchestration. The system uses pre-computed regulatory networks derived from 500,000+ single-cell RNA-seq profiles across 10 cell types from the CellxGene Data Portal, processed through the ARACNe algorithm (a statistical mutual information method for network inference). Four specialized domain agents—cancer biology, drug discovery, clinical relevance, and systems biology—execute in parallel using local language model inference (Ollama/llama3.1:8b) to generate scientific insights with rationales, falling back to rule-based heuristics for reliability. Automated perturbation analysis simulates upstream regulator inhibition for all regulators and ranks therapeutic targets using network centrality metrics (PageRank, degree centrality, out-degree centrality) from computational biology literature. We validated the framework on a colorectal cancer biomarker case study, analyzing five genes (MYC, CTNNB1, CCND1, TP53, KRAS) with focused network and complete perturbation analysis (99 regulators total) completing in 6.15 seconds (rule-based mode) or ~20 seconds (LLM mode with 4 agents × 5 genes); comprehensive analysis including external pathway enrichment (Reactome API) achieves 100% alignment with published literature. Perturbation analysis of TP53 identified seven therapeutic targets ranked by PageRank, successfully predicting experimentally validated Hippo pathway regulators (WWTR1, YAP1) while generating novel hypotheses for experimental prioritization. A conversational interface via Model Context Protocol integration enables natural language access through Claude Desktop, making sophisticated regulatory network analysis accessible without programming expertise. RegNetAgents transforms multi-hour manual workflows into few-second automated analysis while maintaining biological accuracy, providing a scalable framework for gene regulatory research.
```

**Keywords**:
```
gene regulatory networks, multi-agent systems, workflow orchestration, biomarker discovery, therapeutic target identification, network centrality, PageRank, LangGraph, Model Context Protocol, large language models
```

#### Part 5: Author Information

**For each author**:
- First name
- Last name
- Email address
- Institution/Affiliation
- ORCID (optional - get free ID at https://orcid.org)
- Check box: ☑ Corresponding author (if you're the contact person)

**Author contributions** (example):
```
[Your Name]: Conceptualization, methodology, software, validation, writing
```

#### Part 6: Conflict of Interest
- Select: **"No conflicts of interest to declare"**

Or if applicable:
- Describe any conflicts (employment, funding, etc.)

#### Part 7: Ethics
- Question: "Does this study involve human subjects?"
  - Select: **NO**

- Question: "Does this study involve animals?"
  - Select: **NO**

- Question: "Are all data publicly available?"
  - Select: **YES** (CellxGene data, Reactome pathways)

#### Part 8: Funding
- If you have funding, list it here
- If not: "This research received no specific grant funding"

#### Part 9: Data Availability Statement

```
Regulatory network data are derived from publicly available single-cell
RNA-seq datasets in CellxGene (https://cellxgene.cziscience.com/)
processed through the ARACNe algorithm. Pathway annotations use the
Reactome Pathway Database (https://reactome.org). Analysis code is
available from the corresponding author upon reasonable request for
academic research purposes.
```

#### Part 10: License
- Select: **"CC BY 4.0"** (most common, allows others to share/adapt with attribution)
- Or: **"CC BY-NC 4.0"** (non-commercial use only)

---

### Step 4: Review and Submit

1. Review all information carefully
2. Check the PDF preview
3. Read and accept terms and conditions
4. Click **"Submit"**

---

### Step 5: Confirmation

You'll receive:
1. **Immediate**: Confirmation email with submission ID
2. **Within 24-48 hours**: Email saying "under review" by bioRxiv screeners
3. **Within 1-3 days**: Email saying "posted" with your DOI and live link

**DOI format**: Will be something like `doi.org/10.1101/2025.01.20.xxxxxx`

---

## Phase 3: After Posting (~5 minutes)

### Step 1: Update Your Materials

Once you get your DOI and bioRxiv URL:

1. Update `one_page_summary.md`:
   - Line with `[DOI to be added]` → Add your actual DOI

2. Update `email_template.md`:
   - Line with `[bioRxiv URL - add after submission]` → Add your URL

### Step 2: Share (if you want)

**Minimal sharing** (our agreed plan):
- Send 3 emails using the template
- Wait 2 weeks
- Evaluate responses

**Optional sharing**:
- Add to your CV
- List on your lab website
- Mention in your email signature

---

## Common Issues and Solutions

### Issue: PDF upload fails
**Solution**:
- Check file size (max 30 MB)
- Try re-saving as PDF with lower resolution images
- Split very long manuscripts into main + supplementary files

### Issue: "Manuscript doesn't meet criteria"
**Solution**:
- bioRxiv requires new research, not reviews/opinions
- Your work is new research, so this shouldn't happen
- If it does, email bioRxiv support: submit@biorxiv.org

### Issue: Figures aren't showing up
**Solution**:
- Embed figures directly in the PDF before uploading
- Or upload figures as supplementary files

### Issue: Can't find institution in dropdown
**Solution**:
- Type in the search box
- Or select "Other" and type manually

### Issue: Want to update preprint after posting
**Solution**:
- bioRxiv allows revisions
- Log in → "My Submissions" → Upload new version
- Both versions remain available (version tracking)

---

## Timeline Summary

| Day | Activity | Time |
|-----|----------|------|
| **Day 0** | ✅ Prep work (figures, tables, manuscript) | COMPLETE |
| **Today** | Convert Word to PDF, review | 5 minutes |
| **Today** | Submit to bioRxiv | 15-20 minutes |
| **Day 1-3** | bioRxiv screening | (waiting) |
| **Day 3** | Posted online, get DOI | - |
| **Day 3** | Send 3 emails to researchers | 30 minutes |
| **Day 17** | Check responses, decide next steps | 30 minutes |

**Total remaining time**: ~1 hour to submit + 1 hour follow-up over 2 weeks

---

## Checklist

**Already Complete** ✅:
- [x] Edited author names and affiliations in draft
- [x] Generated all 3 figures (PDF and PNG versions)
- [x] Generated all 2 tables (CSV and TXT versions)
- [x] Created Word document with updated author info

**Still To Do**:
- [ ] Convert Word document to PDF (Step 1, Phase 1)
- [ ] Review PDF for errors/sensitive info
- [ ] Create bioRxiv account
- [ ] Submit manuscript and figures (Phase 2)
- [ ] Update DOI links after posting (Phase 3)

---

## Support

**bioRxiv Help**:
- FAQ: https://www.biorxiv.org/about/FAQ
- Email: submit@biorxiv.org
- Response time: 1-2 business days

**Technical Issues**:
- Pandoc: https://pandoc.org/help.html
- Python/matplotlib: Check installation with `python --version` and `pip list`

---

## After Submission

### What bioRxiv Does:
✅ Screens for basic scientific merit (not peer review)
✅ Checks for plagiarism
✅ Assigns DOI
✅ Makes preprint publicly available
✅ Indexes in Google Scholar
✅ Sends email alerts to subscribers in your field

### What You Can Do Next:
✅ Add preprint to your CV
✅ Send to 3 researchers (our minimal plan)
✅ Submit to peer-reviewed journals (preprints don't prevent this)
✅ Post revisions as you get feedback
✅ Track views and downloads in your bioRxiv account

---

## Remember

**This is low risk, low commitment:**
- ✅ FREE to submit
- ✅ Preprints are common and accepted
- ✅ You can update it later
- ✅ It doesn't prevent journal publication
- ✅ You can let it sit quietly or promote it actively (your choice)

**You're doing great!** Take it one step at a time.
