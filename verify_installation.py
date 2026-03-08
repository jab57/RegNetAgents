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
import subprocess
import pickle
import time
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


def check_git_lfs():
    """Check if Git LFS is installed and files are downloaded."""
    print("\nChecking Git LFS...")

    # Check if git lfs command exists
    try:
        result = subprocess.run(['git', 'lfs', 'version'],
                              capture_output=True, text=True, timeout=5)
        if result.returncode == 0:
            version = result.stdout.strip().split('\n')[0]
            print(f"  Git LFS installed: {version}")
        else:
            print(f"{CROSS} Git LFS not installed")
            print("  Install from: https://git-lfs.github.com")
            print("  Or run: git lfs install")
            return False
    except (subprocess.TimeoutExpired, FileNotFoundError):
        print(f"{CROSS} Git LFS not found")
        print("  Install from: https://git-lfs.github.com")
        return False

    # Check if network files are actually downloaded (not just pointers)
    sample_file = Path("models/networks/epithelial_cell/network_index.pkl")

    if not sample_file.exists():
        print(f"{WARN} Sample network file not found")
        print("  Network files may not be present yet")
        return True  # Don't fail - network check will catch this

    # Read first bytes to check if it's a real file or a Git LFS pointer
    try:
        with open(sample_file, 'rb') as f:
            first_bytes = f.read(20)

        # Check if it's a text pointer (Git LFS not pulled)
        if first_bytes.startswith(b'version'):
            print(f"{CROSS} Network files are Git LFS pointers (not downloaded)")
            print("  Fix this by running:")
            print("    git lfs pull")
            print("  This will download ~1.2 GB of network data")
            return False
        # Check if it's a pickle file (correct)
        elif first_bytes.startswith(b'\x80\x03') or first_bytes.startswith(b'\x80\x04'):
            print(f"{CHECK} Network files downloaded correctly")
            return True
        else:
            print(f"{WARN} Network file format unexpected")
            return True  # Let network loading test catch any issues

    except Exception as e:
        print(f"{WARN} Could not verify network file: {e}")
        return True  # Don't fail on read errors


def check_python_version():
    """Verify Python version is 3.10+."""
    print("Checking Python version...")
    version = sys.version_info
    if version.major >= 3 and version.minor >= 10:
        print(f"{CHECK} Python {version.major}.{version.minor}.{version.micro}")
        return True
    else:
        print(f"{CROSS} Python {version.major}.{version.minor}.{version.micro}")
        print("  Required: Python 3.10 or higher")
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
    """Check if network cache files exist and can be loaded."""
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
    total_size = 0
    for cell_type in cell_types:
        cache_file = models_dir / cell_type / "network_index.pkl"
        if cache_file.exists():
            file_size = cache_file.stat().st_size
            # Check if file is reasonable size (> 10KB, not a tiny LFS pointer)
            if file_size > 10_000:  # 10 KB minimum (LFS pointers are ~130 bytes)
                ready_count += 1
                total_size += file_size
            elif file_size < 1000:  # Likely a Git LFS pointer
                print(f"{WARN} {cell_type} network file is too small ({file_size} bytes)")
                print("  This looks like a Git LFS pointer. Run: git lfs pull")

    print(f"  Found {ready_count}/{len(cell_types)} cell type networks ({total_size / 1_000_000:.1f} MB)")

    if ready_count == 0:
        print("[FAIL] No valid network data files found")
        print("  Run: python scripts/build_network_cache.py --all")
        return False

    # Try to load one network to verify it's not corrupted
    print("  Testing network loading...")
    test_file = models_dir / 'epithelial_cell' / 'network_index.pkl'

    if test_file.exists() and test_file.stat().st_size > 10_000:
        try:
            start_time = time.time()
            with open(test_file, 'rb') as f:
                network_data = pickle.load(f)
            load_time = time.time() - start_time

            # Verify it's a dict with expected structure
            if isinstance(network_data, dict) and 'all_genes' in network_data:
                nodes = network_data.get('num_genes', len(network_data.get('all_genes', [])))
                edges = network_data.get('num_edges', 0)
                print(f"{CHECK} Successfully loaded epithelial_cell network ({nodes:,} genes, {edges:,} edges, {load_time:.2f}s)")
            else:
                print(f"{WARN} Network data structure unexpected")
                return True  # Don't fail, but warn

        except Exception as e:
            print(f"{CROSS} Failed to load network: {e}")
            print("  Network files may be corrupted. Try: git lfs pull")
            return False
    else:
        print(f"{WARN} Could not test network loading (epithelial_cell not available)")

    if ready_count == len(cell_types):
        print("[OK] All network data files present and valid")
        return True
    elif ready_count > 0:
        print(f"[WARN] Partial network data ({ready_count}/{len(cell_types)} cell types)")
        print("  Some analyses may not work. Run: git lfs pull")
        return True
    else:
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
        print("[WARN] Ollama Python package not installed (optional)")
        print("  Install with: pip install ollama")
        print("  (Optional - system will use rule-based fallback)")
        return True  # Optional — not required for core functionality


def check_core_modules():
    """Try importing and instantiating core RegNetAgents modules."""
    print("\nChecking core RegNetAgents modules...")

    # Test GeneIDMapper instantiation
    try:
        from regnetagents import GeneIDMapper
        mapper = GeneIDMapper()
        print("[OK] GeneIDMapper (instantiated successfully)")
    except ImportError as e:
        print(f"[FAIL] Failed to import GeneIDMapper: {e}")
        return False
    except Exception as e:
        print(f"[FAIL] Failed to instantiate GeneIDMapper: {e}")
        return False

    # Test CompleteGeneService instantiation
    try:
        from regnetagents import CompleteGeneService
        service = CompleteGeneService()
        print("[OK] CompleteGeneService (instantiated successfully)")
    except ImportError as e:
        print(f"[FAIL] Failed to import CompleteGeneService: {e}")
        return False
    except Exception as e:
        print(f"[WARN] CompleteGeneService import OK but instantiation failed: {e}")
        print("  This is OK if gene annotation file is missing")
        # Don't fail - service can work without annotations

    # Test workflow module import (don't instantiate - requires network data)
    try:
        from regnetagents_langgraph_workflow import RegNetAgentsWorkflow
        print("[OK] LangGraph workflow module (import successful)")
    except ImportError as e:
        print(f"[FAIL] Failed to import workflow module: {e}")
        return False

    # Test MCP server module import
    try:
        from regnetagents_langgraph_mcp_server import server
        print("[OK] MCP server module (import successful)")
    except ImportError as e:
        print(f"[FAIL] Failed to import MCP server: {e}")
        return False

    return True


def check_cache_directory():
    """Check if cache directory exists and is writable."""
    print("\nChecking cache directory...")

    cache_dir = Path("cache")

    # Check if directory exists
    if not cache_dir.exists():
        try:
            cache_dir.mkdir(parents=True, exist_ok=True)
            print(f"{CHECK} Cache directory created: {cache_dir}")
            return True
        except Exception as e:
            print(f"{CROSS} Failed to create cache directory: {e}")
            return False

    # Check if directory is writable
    try:
        test_file = cache_dir / ".write_test"
        test_file.write_text("test")
        test_file.unlink()
        print(f"{CHECK} Cache directory writable: {cache_dir}")
        return True
    except Exception as e:
        print(f"{CROSS} Cache directory not writable: {e}")
        return False


def main():
    """Run all verification checks."""
    print_header("RegNetAgents Installation Verification")

    start_time = time.time()

    checks = [
        ("Python Version", check_python_version()),
        ("Required Packages", check_required_packages()),
        ("Git LFS", check_git_lfs()),
        ("Network Data", check_network_data()),
        ("Cache Directory", check_cache_directory()),
        ("Ollama (optional)", check_ollama()),
        ("Core Modules", check_core_modules())
    ]

    total_time = time.time() - start_time

    print_header("Verification Summary")

    passed = sum(1 for _, result in checks if result)
    total = len(checks)
    required = total - 1  # Ollama is optional

    print(f"\nTotal checks: {total} ({required} required, 1 optional)")
    print(f"Time taken: {total_time:.1f} seconds\n")

    for name, result in checks:
        if result:
            print(f"  {CHECK} {name}")
        else:
            print(f"  {CROSS} {name}")

    print(f"\nResult: {passed}/{total} checks passed")
    print("=" * 60)

    if passed == total:
        print(f"\n{CHECK} SUCCESS! RegNetAgents is ready to use.\n")
        print("Next steps:")
        print("  1. Configure Claude Desktop (see README.md)")
        print("  2. Restart Claude Desktop completely")
        print("  3. Test with: 'Analyze TP53 in epithelial cells'\n")
        return 0
    elif passed >= required:  # All required checks passed (Ollama optional)
        print(f"\n{CHECK} READY TO USE (optional features disabled)\n")
        print("RegNetAgents will work with rule-based mode.")
        print("\nOptional enhancements:")
        print("  • Install Ollama for LLM-powered insights")
        print("    → https://ollama.com/download")
        print("    → ollama pull llama3.1:8b\n")
        return 0
    else:
        print(f"\n{CROSS} INSTALLATION INCOMPLETE\n")
        print("Please fix the failed checks above.")
        print("\nCommon fixes:")
        print("  • Missing packages: pip install -r requirements.txt")
        print("  • Git LFS not pulled: git lfs pull")
        print("  • Network data missing: git lfs install && git lfs pull\n")
        return 1


if __name__ == "__main__":
    sys.exit(main())
