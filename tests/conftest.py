"""Shared pytest configuration for RegNetAgents test suite.

Adds the project root to sys.path so all test files can import
regnetagents modules without per-file path manipulation.
"""

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
