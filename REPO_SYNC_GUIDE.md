# Repository Sync Guide - RegNetAgents

**Reference for syncing backup and public repositories.**

---

## Two Repository Setup

**Backup Repository (RegNetAgents-Backup):**
- Remote name: `backup`
- URL: `https://github.com/jab57/RegNetAgents-Backup.git`
- Purpose: Complete backup of ALL local content
- Visibility: Private
- What goes here: EVERYTHING (code, manuscript, planning docs, all files)

**Public Repository (RegNetAgents):**
- Remote name: `origin`
- URL: `https://github.com/jab57/RegNetAgents.git`
- Purpose: Publication-ready content for bioRxiv/sharing
- Visibility: Public
- What goes here: Only safe/public files (code, manuscript, docs - NO planning)

---

## For Claude AI: Sync Execution Protocol

**When user asks to sync/commit/update repos, follow these exact steps:**

### Step 1: Check Status
```bash
git status
```

### Step 2: Classify Files
Check modified files against the **File Classification** section (lines 331-396):
- **ALL files public-safe?** → Proceed to Step 3a
- **ANY file private?** → Proceed to Step 3b

### Step 3a: Public-Safe Files Only

```bash
# Commit changes
git add <files>
git commit -m "Clear description"

# Push to backup (ALWAYS)
git push backup main

# Try pushing to public
git push origin main
```

**If push succeeds:** ✅ Done!

**If push fails with "non-fast-forward" or "branches have diverged":**
→ Go to **Scenario: Diverged Branches** below

### Step 3b: Private Files Included

```bash
# Commit changes
git add <files>
git commit -m "Clear description"

# Push to backup ONLY
git push backup main

# DO NOT push to origin
```

Output to user: "Pushed to backup only (contains private files: <list files>)"

---

### Scenario: Diverged Branches

**Detection:** `git push origin main` fails with:
```
! [rejected]        main -> main (non-fast-forward)
error: failed to push some refs
```

**Solution - Cherry-pick the commit:**

```bash
# Step 1: Get the commit hash of what you just committed
COMMIT_HASH=$(git rev-parse HEAD)

# Step 2: Checkout origin/main
git checkout origin/main

# Step 3: Cherry-pick the commit
git cherry-pick $COMMIT_HASH

# Step 4: Push to origin
git push origin HEAD:main

# Step 5: Return to main branch
git checkout main
```

**If cherry-pick succeeds:** ✅ Done! Commit is now in both repos.

**If cherry-pick fails with conflicts:**
```
error: could not apply <hash>... <message>
hint: after resolving the conflicts, mark the corrected paths
hint: with 'git add <paths>' or 'git rm <paths>'
```

→ **STOP and ask user:** "Cherry-pick has conflicts. How would you like to proceed?"
→ **DO NOT attempt to resolve automatically**

---

### Scenario: Mixed Public/Private Files

**Detection:** User modified BOTH public-safe AND private files in same working directory

**Solution - Create TWO commits:**

```bash
# Step 1: Stage ONLY public-safe files
git add <public-file-1> <public-file-2>

# Step 2: Commit public files
git commit -m "Description of public changes"

# Step 3: Push to both repos
git push backup main
git push origin main  # (handle divergence if needed - see above)

# Step 4: Stage ONLY private files
git add <private-file-1> <private-file-2>

# Step 5: Commit private files
git commit -m "Description of private changes"

# Step 6: Push to backup ONLY
git push backup main
```

---

### Quick Decision Matrix

| Situation | Files Changed | Action |
|-----------|---------------|--------|
| Normal sync | All public-safe | Push to `backup` AND `origin` |
| Private content | Any private files | Push to `backup` ONLY |
| Push rejected | Diverged branches | Cherry-pick to origin (see above) |
| Mixed changes | Public + private | Split into TWO commits (see above) |
| Cherry-pick conflict | Merge conflict | ASK USER (don't force) |
| Unsure about file | Can't classify | ASK USER or default to backup-only (safer) |

---

### Error Handling

**If unsure which files are public vs private:**
1. Check the **File Classification** section (lines 331-396)
2. If still unsure, ask user: "Is <filename> public-safe or private?"
3. Default to backup-only (safer)

**If git command fails:**
1. Show user the exact error
2. Explain what it means
3. Ask how to proceed (don't guess)

**If push to backup fails:**
1. STOP immediately
2. Alert user - backup is critical
3. Don't proceed with origin push

---

## How It Works (Plain English)

**The Goal:** Keep all your work backed up in one repo, while only sharing safe/public content in another repo.

**The Process:**

1. **You work normally** - Edit any files (code, manuscript, planning docs, anything)

2. **When ready to save** - Tell Claude to commit and backup

3. **Claude checks what changed** - Looks at which files you modified
   ```bash
   $ git status

   # Example output:
   modified:   regnetagents_langgraph_workflow.py
   modified:   biorxiv/preprint_draft.md
   modified:   biorxiv/submission_guide.md
   ```

4. **Claude commits the changes**
   ```bash
   $ git add -A
   $ git commit -m "Fix workflow bug and update docs"

   # Creates a commit with all changes
   ```

5. **Everything goes to backup** - ALL changes get pushed to the private backup repo (RegNetAgents-Backup)
   ```bash
   $ git push backup main

   # ✅ Complete backup of your work
   ```
   - This happens EVERY time, no exceptions

6. **Claude decides about public** - Checks if changed files are safe to share publicly:

   **If ALL files are public-safe** (code, manuscript, docs, figures):
   ```bash
   $ git push origin main

   # ✅ Public repo also gets updated
   ```

   **If ANY file is private** (planning docs, assessments, marketing):
   ```bash
   # ⏸️ Skip pushing to public
   # Only backup has these changes 🔒

   Output: "NOT pushing to public (contains private files)"
   ```

7. **Result:**
   - Backup repo = Complete backup of everything
   - Public repo = Only publication-ready content
   - Private files never leak to public ✅

**Special Situation - When Public is Behind:**

If previous commits had private files, the public repo falls behind. When you later want to update public with safe changes, Claude uses "cherry-picking":

```bash
# Step 1: Check which commits are ahead
$ git log origin/main..main --oneline

# Example output:
5b737ca Enhance REPO_SYNC_GUIDE (safe) ✅
4ea2ed7 Add REPO_SYNC_GUIDE (safe) ✅
0ae4cc5 Add planning documents (PRIVATE) ❌

# Step 2: Cherry-pick only safe commits
$ git cherry-pick 4ea2ed7
$ git cherry-pick 5b737ca

# Step 3: Push to public
$ git push origin main

# Result:
# ✅ Backup has: 0ae4cc5, 4ea2ed7, 5b737ca (everything)
# ✅ Public has: 4ea2ed7, 5b737ca (safe commits only)
# 🔒 Commit 0ae4cc5 never reaches public
```

**You never need to think about this** - just tell Claude to update repos, and the guide ensures everything goes to the right place.

---

## Workflow Schematic

```
┌─────────────────────────────────────────────────────────────┐
│  LOCAL REPOSITORY (C:\Dev\RegNetAgents)                     │
│  Branch: main                                                │
└─────────────────────────────────────────────────────────────┘
                            │
                    User makes changes
                    (edit any files)
                            │
                            ↓
┌─────────────────────────────────────────────────────────────┐
│  STEP 1: Check what changed                                 │
│  $ git status                                                │
│                                                               │
│  Example output:                                             │
│  modified:   regnetagents_langgraph_workflow.py  ← Code     │
│  modified:   biorxiv/preprint_draft.md           ← Manuscript│
│  modified:   biorxiv/submission_guide.md         ← Planning │
└─────────────────────────────────────────────────────────────┘
                            │
                            ↓
┌─────────────────────────────────────────────────────────────┐
│  STEP 2: Commit changes                                     │
│  $ git add -A                                                │
│  $ git commit -m "Fix bug and update docs"                  │
└─────────────────────────────────────────────────────────────┘
                            │
                            ↓
┌─────────────────────────────────────────────────────────────┐
│  STEP 3: ALWAYS push to backup (no exceptions)              │
│  $ git push backup main                                      │
│                                                               │
│  ✅ Everything backed up to RegNetAgents-Backup             │
└─────────────────────────────────────────────────────────────┘
                            │
                            ↓
              ┌─────────────┴─────────────┐
              │  STEP 4: Classify files   │
              └─────────────┬─────────────┘
                            │
              ┌─────────────┴─────────────┐
              │                           │
              ↓                           ↓
    ┌─────────────────┐         ┌─────────────────┐
    │ ALL files are   │         │ ANY file is     │
    │ PUBLIC-SAFE?    │         │ PRIVATE?        │
    │                 │         │                 │
    │ ✅ Code         │         │ ❌ Planning docs│
    │ ✅ Manuscript   │         │ ❌ Assessments  │
    │ ✅ Docs         │         │ ❌ Marketing    │
    │ ✅ Figures      │         │ ❌ Strategy     │
    └────────┬────────┘         └────────┬────────┘
             │                           │
             ↓                           ↓
    ┌─────────────────┐         ┌─────────────────┐
    │ PUSH TO PUBLIC  │         │ SKIP PUBLIC     │
    │ $ git push      │         │ (backup only)   │
    │   origin main   │         │                 │
    └────────┬────────┘         └────────┬────────┘
             │                           │
             ↓                           ↓
┌────────────────────────────┐  ┌────────────────────────────┐
│  BACKUP REPO (Private)     │  │  BACKUP REPO (Private)     │
│  ✅ Has everything         │  │  ✅ Has everything         │
│                            │  │                            │
│  PUBLIC REPO (Public)      │  │  PUBLIC REPO (Public)      │
│  ✅ Also updated           │  │  ⏸️  Not updated          │
│  📢 Visible to all         │  │  🔒 Private files safe    │
└────────────────────────────┘  └────────────────────────────┘


SPECIAL CASE: Cherry-Picking (when public is behind)
═══════════════════════════════════════════════════════════════

If public repo is behind due to previous private commits:

Local commits:
├─ Commit A: Add planning docs (PRIVATE) ❌
├─ Commit B: Fix workflow bug (PUBLIC) ✅
└─ Commit C: Update manuscript (PUBLIC) ✅

Solution: Cherry-pick only public commits
$ git log origin/main..main --oneline    ← See all commits
$ git cherry-pick <commit-B-hash>        ← Pick safe ones
$ git cherry-pick <commit-C-hash>
$ git push origin main                   ← Update public

Result:
- Backup has: A, B, C (everything)
- Public has: B, C only (safe commits)
- Commit A never reaches public ✅
```

---

## Sync Workflow

### Step 1: Check what changed
```bash
git status
```

Review the list of changed files.

### Step 2: Commit changes
```bash
git add -A
git commit -m "Clear description of changes"
```

### Step 3: Always push to backup
```bash
git push backup main
```
**This happens EVERY time, no exceptions.**

### Step 4: Decide about public repo

**IF changes include ONLY these types of files:**
- Code files (*.py)
- README.md, requirements.txt
- Documentation (docs/*.md, except docs/archive/)
- Manuscript (biorxiv/preprint_draft.md)
- Figures (biorxiv/*.png, *.pdf)
- Supplementary materials (biorxiv/supplementary_material.md, SUPPLEMENTARY_MATERIALS.md)
- Tests (tests/*.py)
- Scripts (scripts/*.py)

**THEN also push to public:**
```bash
git push origin main
```

**IF changes include ANY of these files:**
- biorxiv/README.md (internal submission guide)
- biorxiv/submission_guide.md
- biorxiv/journal_submission_plan_REVISED.md
- biorxiv/email_template.md
- biorxiv/PROOFREAD_REPORT.md
- biorxiv/SCIENTIFIC_INTEGRITY_VERIFICATION.md
- biorxiv/timing_verification_protocol.md
- biorxiv/week2_execution_plan.md
- biorxiv/one_page_summary.md
- docs/archive/* (any files)
- Any marketing plans
- Any strategy/planning documents

**THEN do NOT push to public** (backup only).

---

## ⚠️ IMPORTANT: Handling Mixed Changes

**What if you have BOTH public-safe AND private files changed?**

❌ **DO NOT commit them together in one commit!**

✅ **Instead, create TWO separate commits:**

### Step-by-Step for Mixed Changes:

1. **Check what changed:**
   ```bash
   git status
   ```

2. **Identify which files are public vs private** using the File Classification section below

3. **Create TWO commits:**

   **First commit (Public-safe files only):**
   ```bash
   # Stage ONLY public-safe files
   git add README.md biorxiv/preprint_draft.md

   # Commit with clear message
   git commit -m "Fix terminology in manuscript and README"

   # Push to BOTH repos
   git push backup main
   git push origin main
   ```

   **Second commit (Private files only):**
   ```bash
   # Stage ONLY private files
   git add biorxiv/PROOFREAD_REPORT.md biorxiv/one_page_summary.md

   # Commit with clear message
   git commit -m "Update internal review documents"

   # Push to BACKUP ONLY
   git push backup main
   # (Do NOT push to origin)
   ```

4. **Result:**
   - ✅ Public repo gets manuscript updates
   - ✅ Backup repo gets everything
   - 🔒 Private files stay private

### What if you already committed them together?

If you accidentally created one mixed commit, you need to split it:

```bash
# Reset the commit but keep changes
git reset --soft HEAD~1

# Now follow steps 3-4 above to create two separate commits
```

**Key Rule:** One commit = one type (all public OR all private, never mixed)

---

## File Classification

### ✅ PUBLIC-SAFE (can go to both repos)

**Code:**
- `*.py` (all Python files)
- `requirements.txt`
- `.env.example`

**Documentation:**
- `README.md`
- `docs/DATA_SOURCES.md`
- `docs/REGNETAGENTS_MCP_SETUP.md`
- `docs/ADDING_NEW_CELL_TYPES.md`
- `docs/CLAUDE_DESKTOP_USAGE.md`
- `docs/END_TO_END_DATA_PIPELINE.md`
- `docs/GENE_MAPPING_ARCHITECTURE.md`
- `docs/REACTOME_INTEGRATION_SUMMARY.md`
- `docs/REGNETAGENTS_Analysis_Pipeline.md`
- `docs/REGNETAGENTS_CONFERENCE_POSTER.md`

**Manuscript & Publication:**
- `biorxiv/preprint_draft.md`
- `biorxiv/supplementary_material.md`
- `biorxiv/SUPPLEMENTARY_MATERIALS.md`
- `biorxiv/*.png` (all figure images)
- `biorxiv/*.pdf` (all figure PDFs)
- `biorxiv/table*.csv` (data tables)
- `biorxiv/table*.txt` (data tables)
- `biorxiv/CONVERSION_README.md`
- `biorxiv/create_figure1.py`
- `biorxiv/convert_manuscript.py`
- `biorxiv/regnetagents_preprint.docx`
- `biorxiv/regnetagents_preprint.pdf`

**Data:**
- `models/networks/*` (network cache files)
- `cache/gene_id_cache.pkl`

**Tests:**
- `tests/*.py`
- `scripts/*.py`

### ❌ PRIVATE (backup only, NOT public)

**Planning Documents:**
- `biorxiv/README.md` (internal submission guide)
- `biorxiv/submission_guide.md`
- `biorxiv/journal_submission_plan_REVISED.md`
- `biorxiv/email_template.md`
- `biorxiv/week2_execution_plan.md`
- `biorxiv/one_page_summary.md`

**Internal Reviews:**
- `biorxiv/PROOFREAD_REPORT.md`
- `biorxiv/SCIENTIFIC_INTEGRITY_VERIFICATION.md`
- `biorxiv/timing_verification_protocol.md`

**Internal Assessments:**
- `docs/archive/*` (all files in archive)

**Strategy/Marketing:**
- Any marketing plans
- Any strategy documents
- Any internal business planning

---

## Decision Tree

```
Changed files?
    │
    ├─ ONLY public-safe files?
    │   └─ YES → Push to backup AND origin
    │
    └─ Includes ANY private files?
        └─ YES → Push to backup ONLY (not origin)
```

---

## Commands Reference

**Check status:**
```bash
git status
```

**Commit changes:**
```bash
git add -A
git commit -m "Description of changes"
```

**Push to backup (always):**
```bash
git push backup main
```

**Push to public (conditional):**
```bash
git push origin main
```

**Push to both (when safe):**
```bash
git push backup main && git push origin main
```

---

## Current State Check

**See what's in each repo:**
```bash
# Public repo commits
git log origin/main --oneline -5

# Backup repo commits
git log backup/main --oneline -5

# Commits in backup but not public
git log origin/main..backup/main --oneline
```

---

## Notes

- Documentation already points to public repo (github.com/jab57/RegNetAgents)
- Backup repo contains complete backup (cannot rely on OneDrive/external drives)
- Some planning docs may contain marketing/strategy info (keep private)
- When in doubt, push to backup only (safer to be conservative)

---

## Summary

1. **Check** what files changed
2. **Commit** with clear message
3. **Always** push to backup
4. **Conditionally** push to public (only if all files are public-safe)
5. **When uncertain** → backup only
