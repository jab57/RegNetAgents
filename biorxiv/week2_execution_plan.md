# Week 2: Framework Value Demonstration

## Goal
Create Figure 4 showing RegNetAgents is **orders of magnitude faster** than manual workflows and **LLM adds scientific context**

## Timeline: 4-6 days

---

## Tasks

### 1. Create Figure 4 (Days 1-3)
**File:** `scripts/generate_figures.py` (add after line 459)
**Output:** `biorxiv/figure4_framework_value.png` + `.pdf`

#### Figure 4 Layout (2×2):

### Figure 4 Layout: 2×2 Grid

```
┌─────────────────────┬─────────────────────┐
│  A) Manual Workflow │ B) RegNetAgents     │
│     Flowchart       │    Flowchart        │
│     2-4 hours       │    0.6-15 sec       │
├─────────────────────┼─────────────────────┤
│  C) Time Comparison │ D) LLM Intelligence │
│     Bar Chart       │    Side-by-Side     │
│  orders of magnitude faster │   Scores vs Context │
└─────────────────────┴─────────────────────┘
```

- **Panel A:** Manual workflow (4 steps, 2-4 hours) - gray boxes
- **Panel B:** RegNetAgents workflow (1 query, 15 sec) - green boxes
- **Panel C:** Bar chart comparing times (use data from `results/comprehensive_timing_results.json`)
- **Panel D:** LLM vs rule-based (side-by-side text boxes)

**Data:** All measurements already in `comprehensive_timing_results.json`

### 2. Add Manuscript Sections (Days 4-5)
**File:** `biorxiv/preprint_draft.md`

**Insert after line 265:**
- Section: "Workflow Automation and Performance" (~350 words)
- Section: "LLM-Powered Domain Analysis" (~450 words)

**Insert after line 568:**
- Figure 4 legend (~200 words)

**Update line ~16:**
- Add: "representing orders of magnitude speedup over manual workflows (Figure 4)"

### 3. Validation (Day 6)
- Check Figure 4 renders correctly (300 DPI)
- Verify timing data matches JSON file
- Confirm speedup is orders of magnitude (manual baseline: 2.5 hours vs RegNetAgents: 0.6-15 sec)
- Proofread new manuscript sections

---

## Files to Modify
1. `scripts/generate_figures.py` - Add `create_figure4()` after line 459
2. `biorxiv/preprint_draft.md` - Lines 16, 265, 568
3. `scripts/generate_conference_poster.py` - Add Figure 4 to poster layout

## Success: Week 2 Complete
✅ Figure 4 created (300 DPI, PNG + PDF)
✅ ~800 new words added to manuscript
✅ Framework value demonstrated (speed + intelligence)
