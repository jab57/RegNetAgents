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
- ✅ Your affiliation shows: "Bird AI Solutions"

**Note**: The Word document doesn't include embedded figures. You'll upload figures separately to bioRxiv as supplementary files (instructions in Phase 2).

---

### Step 2: Make GitHub Repository Public (~2 minutes)

Before submitting to bioRxiv, make your GitHub repository public so reviewers and readers can access the code.

**Steps:**

1. Go to https://github.com/jab57/RegNetAgents
2. Click **Settings** (top right, gear icon)
3. Scroll down to **Danger Zone** (bottom of page)
4. Click **Change repository visibility**
5. Select **Make public**
6. Type the repository name to confirm: `jab57/RegNetAgents`
7. Click **I understand, make this repository public**

**What gets shared:**
- ✅ All code and documentation
- ✅ Commit history
- ✅ Issues and pull requests
- ❌ NOT shared: Any files in `.gitignore` (local configs, cache files)

**Important:** Review your repository one last time before making public:
- [ ] No API keys or passwords in code
- [ ] No sensitive data files
- [ ] README.md is clear and informative

**After making public:**
- Add the repository URL to your bioRxiv submission metadata
- Include it in the Data Availability Statement: "Code available at https://github.com/jab57/RegNetAgents"

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

#### Part 3b: Upload Supplementary Files

Look for "Supplementary Files" or "Additional Files" section. Upload files in this order:

**1. Supplementary Material (FIRST):**
- File: `supplementary_material.pdf`
- Label: "Supplementary Material"
- Description: "LLM Prompt Templates and Reproducibility Notes"

**2. Figures:**
- `figure1_architecture.pdf` → Label: "Figure 1"
- `figure2_biomarker_panel.pdf` → Label: "Figure 2"
- `figure3_tp53_perturbation.pdf` → Label: "Figure 3"

**Total supplementary files:** 4 (1 supplementary material + 3 figures)

#### Part 4: Manuscript Information

**Title** (copy from your PDF):
```
RegNetAgents: A Multi-Agent AI Framework for Automated Gene Regulatory Network Analysis and Therapeutic Target Prioritization
```

**Running Title**:
```
Multi-Agent Gene Regulatory Network Analysis
```

**Abstract** (copy from the updated preprint draft):
```
Gene regulatory network analysis is essential for understanding disease mechanisms, identifying biomarkers, and prioritizing therapeutic targets. However, traditional workflows require hours of manual effort across multiple databases and tools, with sequential processing limiting scalability. We present RegNetAgents, an LLM-powered multi-agent framework that streamlines gene regulatory analysis through intelligent workflow orchestration and conversational interfaces. The system integrates network modeling, network-based target ranking, pathway enrichment, and multi-domain interpretation (cancer, drug development, clinical relevance, systems biology) into a unified accessible platform. Leveraging pre-computed ARACNe networks from 500,000+ single cells across 10 cell types (GREmLN team), the framework enables rapid hypothesis generation and experimental prioritization. Four specialized domain agents execute in parallel using local language models to generate scientific insights, with rule-based fallback for reliability. Automated perturbation analysis ranks regulators by network centrality and ranks candidate therapeutic targets using network centrality metrics (PageRank, degree centrality). To demonstrate framework capabilities, we analyzed a colorectal cancer biomarker panel (MYC, CTNNB1, CCND1, TP53, KRAS) with complete perturbation analysis (99 regulators) completing in 15-62 seconds depending on mode. Framework demonstration on colorectal cancer biomarkers showed complete concordance with published literature across this limited five-gene sample, and perturbation analysis identified literature-supported TP53 interactors (WWTR1, YAP1 from Hippo pathway) alongside novel testable hypotheses. A conversational interface via Model Context Protocol enables natural language queries through Claude Desktop without programming expertise. RegNetAgents transforms multi-hour manual workflows into second-scale automated analysis, making sophisticated network analysis accessible to experimental biologists and providing a reusable framework for diverse biological questions.
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
Reactome Pathway Database (https://reactome.org). Analysis code and
implementation details are publicly available at
https://github.com/jab57/RegNetAgents under an open-source license.
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
- [x] Moved appendix to supplementary material
- [x] Generated supplementary_material.pdf (263 KB, ~9 pages)
- [x] Updated main manuscript to reference supplementary material
- [x] Regenerated main manuscript PDF (356 KB, ~40 pages)

**Still To Do**:
- [ ] Review main manuscript PDF for errors
- [ ] Review supplementary material PDF for formatting
- [ ] Make GitHub repository public (Step 2, Phase 1)
- [ ] Create bioRxiv account
- [ ] Submit manuscript, supplementary material, and figures (Phase 2)
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
