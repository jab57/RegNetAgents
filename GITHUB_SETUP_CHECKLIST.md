# GitHub Setup Checklist (For Repository Owner)

This is YOUR checklist for preparing and pushing RegNetAgents to a private GitHub repository.

**Estimated time:** 30-45 minutes

---

## Step 1: Clean Up Local Repository (10 minutes)

### 1.1 Update .gitignore

Add these entries to `.gitignore`:

```
# Environment variables (SENSITIVE - never commit!)
.env

# GREmLN source (separate repository)
gremln_source/

# Results and outputs (optional - users generate their own)
results/*.json

# Temporary files
*.tmp
*.log
.pytest_cache/

# OS-specific
.DS_Store
Thumbs.db
```

**Command:**
```bash
# Open .gitignore and add the lines above, or I can do this for you
```

### 1.2 Review Sensitive Files

**Check these files for sensitive information:**
- [ ] `.env` - Should NOT be committed (add to .gitignore)
- [ ] `.env.example` - Should be committed (template only, no real secrets)
- [ ] Any API keys or tokens in code

**Action:** If `.env.example` doesn't exist, create it:
```bash
# .env.example (template for users)
OLLAMA_HOST=http://localhost:11434
OLLAMA_MODEL=llama3.1:8b
```

---

## Step 2: Commit Current Changes (5 minutes)

You have many pending changes. Let's commit them:

```bash
cd C:/Dev/RegNetAgents

# Stage all changes
git add .

# Commit with descriptive message
git commit -m "Refactor: Rename from GREmLN to RegNetAgents

- Rename main modules to regnetagents_*
- Reorganize into regnetagents/ package
- Move tests to tests/ directory
- Update documentation and README
- Add bioRxiv publication materials"

# Check status
git status
```

**Expected result:** Clean working directory (no uncommitted changes)

---

## Step 3: Verify Git LFS (2 minutes)

Git LFS is already tracking your large data files. Verify it's configured:

```bash
# Check LFS is installed
git lfs version

# See what's tracked
git lfs ls-files

# Should show:
# - data/epithelial_cells.h5ad
# - data/human_immune_cells.h5ad
# - models/model.ckpt
```

**If Git LFS is not installed:**
- Windows: `winget install git-lfs` or download from https://git-lfs.com
- After install: `git lfs install`

---

## Step 4: Create Private GitHub Repository (5 minutes)

### 4.1 Create Repository

1. Go to https://github.com
2. Click **"New"** (green button, top right)
3. Repository name: `RegNetAgents`
4. Description: `Multi-agent AI framework for gene regulatory network analysis via Model Context Protocol`
5. **Visibility: Private** ⚠️
6. **DO NOT** check "Initialize with README" (you already have one)
7. Click **"Create repository"**

### 4.2 Copy the Repository URL

GitHub will show you the URL. It looks like:
```
https://github.com/[your-username]/RegNetAgents.git
```

**Keep this handy for the next step.**

---

## Step 5: Push to GitHub (10 minutes)

### 5.1 Add Remote

```bash
cd C:/Dev/RegNetAgents

# Add GitHub as remote (replace [your-username] with your GitHub username)
git remote add origin https://github.com/[your-username]/RegNetAgents.git

# Verify
git remote -v
```

### 5.2 Push Code

```bash
# Push main branch
git push -u origin main

# This will take a few minutes due to LFS files (171MB + 115MB + 33MB)
```

**Expected output:**
- Uploading objects
- Uploading LFS objects (progress bars for large files)
- Success message

### 5.3 Verify on GitHub

1. Go to your GitHub repository URL
2. Check that you see:
   - ✅ Source code files
   - ✅ README.md displays properly
   - ✅ LFS files show file size (not content)
   - ✅ Repository shows "Private" badge

---

## Step 6: Set Up Repository Settings (5 minutes)

### 6.1 Add Repository Topics

On GitHub repository page:
1. Click **"About"** (gear icon)
2. Add topics: `gene-regulatory-networks`, `mcp-server`, `claude-code`, `bioinformatics`, `multi-agent-systems`, `langgraph`
3. Save changes

### 6.2 Configure Default Branch

1. Go to **Settings** → **Branches**
2. Ensure default branch is `main`
3. (Optional) Add branch protection rules later

---

## Step 7: Test Clone (Optional, 5 minutes)

Test that someone can clone your repository:

```bash
# Clone to a different location to test
cd C:/Temp
git clone https://github.com/[your-username]/RegNetAgents.git RegNetAgents-test

# Check it worked
cd RegNetAgents-test
ls
```

**Clean up test clone:**
```bash
cd ..
rm -rf RegNetAgents-test
```

---

## Step 8: Create Release (Optional, Future)

When ready to share more widely:

1. Go to repository → **Releases** → **Create a new release**
2. Tag: `v0.1.0`
3. Title: `Initial Release - RegNetAgents v0.1.0`
4. Description: Brief overview
5. Mark as "pre-release" if still testing
6. Publish

---

## Troubleshooting

### Issue: LFS push fails

**Solution:**
```bash
# Check LFS bandwidth quota
git lfs env

# If quota exceeded, wait or upgrade GitHub plan
```

### Issue: Push rejected (non-fast-forward)

**Solution:**
```bash
# Pull first (if you made changes on GitHub)
git pull --rebase origin main
git push -u origin main
```

### Issue: .env file accidentally committed

**Solution:**
```bash
# Remove from git but keep local file
git rm --cached .env
git commit -m "Remove .env from repository"
git push

# Then update .gitignore to include .env
```

---

## Checklist Summary

Before considering this complete, verify:

- [ ] .gitignore updated with .env, gremln_source/, results/
- [ ] .env.example exists (template without secrets)
- [ ] All changes committed locally
- [ ] Git LFS is installed and tracking large files
- [ ] Private GitHub repository created
- [ ] Code pushed to GitHub successfully
- [ ] LFS files uploaded (check repository size ~320MB)
- [ ] README displays correctly on GitHub
- [ ] Repository is marked "Private"

---

## Next Steps

After completing this checklist:

1. **Users can now clone** - Share repository URL with collaborators
2. **They need INSTALL.md** - Guide them through setup (see INSTALL.md)
3. **Consider CI/CD** - Add GitHub Actions for testing (future)
4. **Track issues** - Use GitHub Issues for bug tracking
5. **Plan release** - When ready to make public, change visibility

---

## Support

**GitHub Help:**
- Git LFS: https://git-lfs.com
- GitHub authentication: https://docs.github.com/en/authentication
- Private repositories: https://docs.github.com/en/repositories

**Need help?** Review git status and error messages carefully, or consult GitHub documentation.
