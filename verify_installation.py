#!/usr/bin/env python3
"""
RegNetAgents Installation Verification Script
==============================================

This script verifies that RegNetAgents is properly installed and configured.
Run this after installation to ensure all components are working correctly.

Usage:
    python verify_installation.py

Author: Jose A. Bird, PhD
License: MIT
"""

import sys
import os
from pathlib import Path

# Use ASCII-compatible symbols for cross-platform compatibility
CHECK = "[OK]"
CROSS = "[FAIL]"
WARN = "[WARN]"


def print_header(text):
    """Print a formatted header."""
    print(f"\n{'=' * 60}")
    print(f"  {text}")
    print(f"{'=' * 60}\n")


def check_python_version():
    """Verify Python version is 3.8+."""
    print("Checking Python version...")
    version = sys.version_info
    if version.major >= 3 and version.minor >= 8:
        print(f"{CHECK} Python {version.major}.{version.minor}.{version.micro}")
        return True
    else:
        print(f"{CROSS} Python {version.major}.{version.minor}.{version.micro}")
        print("  Required: Python 3.8 or higher")
        return False


def check_required_packages():
    """Check if required packages are installed."""
    print("\nChecking required packages...")

    required_packages = [
        ('langgraph', 'LangGraph'),
        ('mcp', 'Model Context Protocol'),
        ('networkx', 'NetworkX'),
        ('pandas', 'Pandas'),
        ('numpy', 'NumPy'),
        ('requests', 'Requests'),
        ('matplotlib', 'Matplotlib'),
        ('seaborn', 'Seaborn'),
        ('ollama', 'Ollama Python SDK'),
        ('dotenv', 'Python-dotenv')
    ]

    all_installed = True
    for package, name in required_packages:
        try:
            if package == 'dotenv':
                __import__('dotenv')
            else:
                __import__(package)
            print(f"[OK] {name}")
        except ImportError:
            print(f"[FAIL] {name} (NOT INSTALLED)")
            all_installed = False

    return all_installed


def check_network_data():
    """Check if network cache files exist."""
    print("\nChecking network data files...")

    cell_types = [
        'epithelial_cell',
        'cd14_monocytes',
        'cd16_monocytes',
        'cd20_b_cells',
        'cd4_t_cells',
        'cd8_t_cells',
        'erythrocytes',
        'nk_cells',
        'nkt_cells',
        'monocyte-derived_dendritic_cells'
    ]

    models_dir = Path("models/networks")

    if not models_dir.exists():
        print(f"[FAIL] Models directory not found: {models_dir}")
        print("  Please ensure network data is available")
        return False

    ready_count = 0
    for cell_type in cell_types:
        cache_file = models_dir / cell_type / "network_index.pkl"
        if cache_file.exists():
            ready_count += 1

    print(f"  Found {ready_count}/{len(cell_types)} cell type networks")

    if ready_count == len(cell_types):
        print("[OK] All network data files present")
        return True
    elif ready_count > 0:
        print(f"[WARN] Partial network data ({ready_count}/{len(cell_types)} cell types)")
        return True
    else:
        print("[FAIL] No network data files found")
        print("  Run: python scripts/build_network_cache.py --all")
        return False


def check_ollama():
    """Check if Ollama is available (optional)."""
    print("\nChecking Ollama (optional, for LLM agents)...")

    try:
        import ollama
        try:
            # Try to list models to verify Ollama is running
            models = ollama.list()
            model_names = [m['name'] for m in models.get('models', [])]

            if 'llama3.1:8b' in model_names:
                print("[OK] Ollama running with llama3.1:8b model")
                return True
            elif model_names:
                print(f"[WARN] Ollama running but llama3.1:8b not found")
                print(f"  Available models: {', '.join(model_names)}")
                print("  Run: ollama pull llama3.1:8b")
                return True
            else:
                print("[WARN] Ollama running but no models installed")
                print("  Run: ollama pull llama3.1:8b")
                return True
        except Exception as e:
            print(f"[WARN] Ollama not running: {e}")
            print("  Install from: https://ollama.com/download")
            print("  (Optional - system will use rule-based fallback)")
            return True  # Still OK, just no LLM
    except ImportError:
        print("[FAIL] Ollama Python package not installed")
        return False


def check_core_modules():
    """Try importing core RegNetAgents modules."""
    print("\nChecking core RegNetAgents modules...")

    try:
        from regnetagents import GeneIDMapper, CompleteGeneService
        print("[OK] RegNetAgents core modules")
    except ImportError as e:
        print(f"[FAIL] Failed to import RegNetAgents modules: {e}")
        return False

    try:
        from regnetagents_langgraph_workflow import RegNetAgentsWorkflow
        print("[OK] LangGraph workflow module")
    except ImportError as e:
        print(f"[FAIL] Failed to import workflow module: {e}")
        return False

    try:
        from regnetagents_langgraph_mcp_server import server
        print("[OK] MCP server module")
    except ImportError as e:
        print(f"[FAIL] Failed to import MCP server: {e}")
        return False

    return True


def main():
    """Run all verification checks."""
    print_header("RegNetAgents Installation Verification")

    checks = [
        ("Python Version", check_python_version()),
        ("Required Packages", check_required_packages()),
        ("Network Data", check_network_data()),
        ("Ollama (optional)", check_ollama()),
        ("Core Modules", check_core_modules())
    ]

    print_header("Verification Summary")

    passed = sum(1 for _, result in checks if result)
    total = len(checks)

    for name, result in checks:
        status = "[OK] PASS" if result else "[FAIL] FAIL"
        print(f"{status:10} | {name}")

    print(f"\n{passed}/{total} checks passed\n")

    if passed == total:
        print("SUCCESS! SUCCESS! RegNetAgents is properly installed.")
        print("\nNext steps:")
        print("1. Configure Claude Desktop (see README.md Step 5)")
        print("2. Restart Claude Desktop")
        print("3. Try: 'Analyze TP53 in epithelial cells'")
        return 0
    elif passed >= total - 1:  # Allow Ollama to be missing
        print("[WARN] PARTIAL SUCCESS. RegNetAgents is mostly configured.")
        print("\nYou can use RegNetAgents, but some features may be limited.")
        print("Check the failed items above for details.")
        return 1
    else:
        print("[FAIL] INSTALLATION INCOMPLETE. Please fix the failed checks above.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
