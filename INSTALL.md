# RegNetAgents Installation Guide

Complete installation instructions for setting up RegNetAgents MCP Server with Claude Desktop.

**Estimated time:** 20-30 minutes (first-time setup)

---

## Prerequisites

Before you begin, ensure you have:

- **Python 3.10 or later** ([Download](https://www.python.org/downloads/))
- **Git** ([Download](https://git-scm.com/downloads))
- **Git LFS** ([Download](https://git-lfs.com)) - Required for large data files
- **Claude Desktop** ([Download](https://claude.ai/download))
- **8GB+ RAM** recommended (for running local LLM)
- **5GB+ free disk space** (for Ollama models)

**Optional but recommended:**
- **Ollama** for LLM-powered domain analysis ([Download](https://ollama.com/download))

---

## Quick Start (TL;DR)

For experienced users:

```bash
# 1. Clone repository
git clone https://github.com/jab57/RegNetAgents.git
cd RegNetAgents

# 2. Install dependencies
pip install -r requirements.txt

# 3. Install Ollama and pull model
ollama pull llama3.1:8b

# 4. Configure Claude Desktop (see Step 4)

# 5. Test
python regnetagents_langgraph_mcp_server.py
```

For detailed instructions, continue reading below.

---

## Step 1: Clone Repository (2 minutes)

### 1.1 Install Git LFS

Git LFS is required to download large data files (network datasets).

**Windows:**
```bash
winget install git-lfs
git lfs install
```

**macOS:**
```bash
brew install git-lfs
git lfs install
```

**Linux:**
```bash
sudo apt-get install git-lfs  # Debian/Ubuntu
# or
sudo yum install git-lfs      # RedHat/CentOS
git lfs install
```

### 1.2 Clone the Repository

```bash
# Clone repository
git clone https://github.com/jab57/RegNetAgents.git

# Navigate to directory
cd RegNetAgents

# Verify LFS files downloaded (should show ~320MB)
git lfs ls-files
```

**Expected LFS files:**
- `data/epithelial_cells.h5ad` (33 MB)
- `data/human_immune_cells.h5ad` (171 MB)
- `models/model.ckpt` (115 MB)

---

## Step 2: Install Python Dependencies (5 minutes)

### 2.1 Create Virtual Environment (Recommended)

```bash
# Create virtual environment
python -m venv regnetagents-env

# Activate virtual environment
# Windows:
regnetagents-env\Scripts\activate

# macOS/Linux:
source regnetagents-env/bin/activate
```

### 2.2 Install Requirements

```bash
# Install all dependencies
pip install -r requirements.txt

# Verify installation
pip list | grep langgraph
pip list | grep mcp
```

**Expected packages:**
- langgraph (workflow orchestration)
- mcp (Model Context Protocol)
- networkx (network analysis)
- pandas, numpy (data processing)
- requests (API calls)
- ollama (LLM inference, optional)

---

## Step 3: Install Ollama (Optional, 10 minutes)

Ollama enables LLM-powered domain analysis with scientific rationales. If you skip this, RegNetAgents will use rule-based analysis (still fully functional).

### 3.1 Install Ollama

**Windows / macOS:**
1. Download from https://ollama.com/download
2. Run installer
3. Ollama runs as background service automatically

**Linux:**
```bash
curl -fsSL https://ollama.com/install.sh | sh
```

### 3.2 Pull Required Model

```bash
# Download llama3.1:8b model (~4.7 GB)
ollama pull llama3.1:8b

# Verify installation
ollama list
```

**Expected output:**
```
NAME              ID              SIZE
llama3.1:8b       a7f6c5...       4.7 GB
```

### 3.3 Test Ollama

```bash
# Test inference
ollama run llama3.1:8b "Hello, are you working?"

# Should respond with a greeting
```

**Note:** First run may take 10-20 seconds to load model into memory.

---

## Step 4: Configure Claude Desktop MCP (5 minutes)

### 4.1 Create Environment File

Create a `.env` file in the RegNetAgents directory:

```bash
# Copy template
cp .env.example .env

# Edit .env (optional - defaults work for most users)
```

**.env contents:**
```
# Ollama configuration
OLLAMA_HOST=http://localhost:11434
OLLAMA_MODEL=llama3.1:8b
```

### 4.2 Configure Claude Desktop

**Location of MCP settings file:**

- **Windows:** `%APPDATA%\Claude\claude_desktop_config.json`
- **macOS:** `~/Library/Application Support/Claude/claude_desktop_config.json`
- **Linux:** `~/.config/Claude/claude_desktop_config.json`

**Edit the file and add RegNetAgents server:**

```json
{
  "mcpServers": {
    "regnetagents": {
      "command": "python",
      "args": [
        "C:/Dev/RegNetAgents/regnetagents_langgraph_mcp_server.py"
      ],
      "env": {
        "PYTHONPATH": "C:/Dev/RegNetAgents"
      }
    }
  }
}
```

**Important:**
- Replace `C:/Dev/RegNetAgents` with your actual installation path
- Use forward slashes `/` even on Windows
- If you have other MCP servers, add regnetagents to the existing list

### 4.3 Restart Claude Desktop

1. Quit Claude Desktop completely
2. Relaunch Claude Desktop
3. Look for "RegNetAgents" in the MCP servers list (bottom right or settings)

---

## Step 5: Verify Installation (5 minutes)

### 5.1 Test MCP Server Standalone

```bash
# Activate virtual environment if not already active
# Windows: regnetagents-env\Scripts\activate
# macOS/Linux: source regnetagents-env/bin/activate

# Run server test
python regnetagents_langgraph_mcp_server.py
```

**Expected output:**
```
MCP server initialized
Registered tools: ['query_gene_network', 'perturbation_analysis', ...]
Server running on stdio...
```

Press `Ctrl+C` to stop.

### 5.2 Test in Claude Desktop

Open Claude Desktop and try these queries:

**Query 1: Simple network lookup**
```
What regulates the TP53 gene in epithelial cells?
```

**Query 2: Perturbation analysis**
```
Analyze therapeutic targets for TP53 in epithelial cells
```

**Query 3: Multi-gene analysis**
```
Compare MYC, TP53, and KRAS across all cell types
```

**Expected behavior:**
- Claude should invoke RegNetAgents tools
- Results appear within seconds (rule-based) or ~15-60 seconds (LLM-powered)
- Output includes network metrics, regulators, targets, and insights

---

## Troubleshooting

### Issue: Git LFS files not downloading

**Symptoms:** Data files show as text pointers instead of binary data

**Solution:**
```bash
git lfs install
git lfs pull
```

### Issue: Import errors (ModuleNotFoundError)

**Symptoms:** `ModuleNotFoundError: No module named 'langgraph'`

**Solution:**
```bash
# Ensure virtual environment is active
pip install -r requirements.txt

# Verify installation
pip list
```

### Issue: Ollama not responding

**Symptoms:** `Connection refused` or timeout errors

**Solution:**
```bash
# Check if Ollama is running
ollama list

# Restart Ollama service
# Windows/macOS: Restart from system tray
# Linux:
sudo systemctl restart ollama
```

### Issue: Claude Desktop doesn't see MCP server

**Symptoms:** RegNetAgents not listed in Claude Desktop

**Solution:**
1. Check `claude_desktop_config.json` syntax (valid JSON)
2. Verify file path is correct (use forward slashes)
3. Restart Claude Desktop completely (quit, not minimize)
4. Check Claude Desktop logs:
   - Windows: `%APPDATA%\Claude\logs\`
   - macOS: `~/Library/Logs/Claude/`

### Issue: LLM analysis not working

**Symptoms:** Results show `llm_powered: false`

**Solution:**
```bash
# Verify Ollama is running
ollama list

# Test model
ollama run llama3.1:8b "test"

# Check .env file has correct settings
cat .env
```

**Note:** If Ollama is unavailable, RegNetAgents automatically falls back to rule-based analysis (still fully functional).

### Issue: Slow performance

**Symptoms:** Queries take >60 seconds

**Possible causes:**
1. **First run:** Network caches loading for first time (~2-3 seconds)
2. **LLM cold start:** Ollama loading model into memory (~10-20 seconds first query)
3. **Comprehensive mode:** Running 4 domain agents in parallel takes 15-60 seconds (expected)

**Solutions:**
- Use focused mode for faster results: "Quick analysis of TP53"
- Increase RAM if system is swapping
- Disable LLM mode if only need network metrics

---

## System Requirements

### Minimum:
- Python 3.10+
- 4GB RAM
- 2GB disk space (without Ollama)
- Network connection (for initial setup, API calls)

### Recommended:
- Python 3.11+
- 8GB+ RAM (for Ollama)
- 10GB disk space (with Ollama models)
- SSD for faster cache loading

### Tested Platforms:
- ✅ Windows 10/11
- ✅ macOS 12+ (Intel & Apple Silicon)
- ✅ Ubuntu 20.04+ / Debian 11+

---

## What's Installed?

After completing installation, you have:

**Core components:**
- RegNetAgents Python package (`regnetagents/`)
- MCP server (`regnetagents_langgraph_mcp_server.py`)
- LangGraph workflow engine (`regnetagents_langgraph_workflow.py`)

**Data:**
- 10 cell-type-specific regulatory networks (models/networks/)
- Pre-computed network caches for fast lookup
- Gene ID mapping tables

**Optional:**
- Ollama with llama3.1:8b model (4.7GB)
- 4 specialized domain analysis agents (cancer, drug, clinical, systems)

---

## Usage Examples

Once installed, try these example queries in Claude Desktop:

### Basic Network Analysis
```
What genes regulate APC in epithelial cells?
```

### Perturbation Analysis
```
Find therapeutic targets for MYC inhibition in epithelial cells
```

### Cross-Cell-Type Comparison
```
How does TP53 behave differently across immune cells vs epithelial cells?
```

### Pathway Enrichment
```
What pathways are enriched for IL6 and its regulators?
```

### Multi-Gene Biomarker Panel
```
Analyze these colorectal cancer biomarkers: MYC, CTNNB1, TP53, KRAS
```

---

## Performance Benchmarks

**Expected execution times:**

| Analysis Type | Rule-Based | LLM-Powered |
|---------------|------------|-------------|
| Single gene (focused) | 0.5-1 sec | 15-20 sec |
| Single gene (comprehensive) | 0.5-1 sec | 15-20 sec |
| Multi-gene (5 genes) | 5-10 sec | 60-90 sec |
| Perturbation analysis | 1-2 sec | 20-30 sec |

**Network cache loading:** <100ms (after first load)

---

## Optional: Running Tests

To verify everything works correctly:

```bash
# Run all tests
pytest tests/

# Run specific test
pytest tests/test_langgraph_workflow.py

# Run with verbose output
pytest -v tests/
```

**Expected result:** All tests pass ✅

---

## Updating

To get the latest version:

```bash
cd RegNetAgents
git pull
pip install -r requirements.txt --upgrade
```

---

## Uninstalling

To completely remove RegNetAgents:

```bash
# 1. Remove from Claude Desktop config
# Edit claude_desktop_config.json and remove "regnetagents" entry

# 2. Delete virtual environment
rm -rf regnetagents-env/

# 3. Delete repository
cd ..
rm -rf RegNetAgents/

# 4. (Optional) Remove Ollama model
ollama rm llama3.1:8b
```

---

## Getting Help

**Documentation:**
- See `README.md` for overview and features
- See `docs/` directory for detailed guides
- See `TROUBLESHOOTING.md` for common issues

**Support:**
- GitHub Issues: [repository URL]/issues
- Email: [your email]

**Contributing:**
- See `CONTRIBUTING.md` for guidelines (if available)

---

## Next Steps

After successful installation:

1. **Read the README** - Understand features and capabilities
2. **Try example queries** - Test with sample analyses
3. **Explore documentation** - See `docs/` for advanced usage
4. **Run your analyses** - Apply to your research questions
5. **Provide feedback** - Report bugs or suggest features

---

**Installation complete!** 🎉

You're ready to use RegNetAgents for gene regulatory network analysis through Claude Desktop.
