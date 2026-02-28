# Contributing to RegNetAgents

Thank you for your interest in contributing to RegNetAgents.

## How to Contribute

### Reporting Issues

- Use [GitHub Issues](https://github.com/jab57/RegNetAgents/issues) to report bugs or request features
- Include your Python version, operating system, and steps to reproduce

### Submitting Changes

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/your-feature`)
3. Install development dependencies: `pip install -e ".[dev]"`
4. Make your changes
5. Run tests: `pytest tests/ --cov=regnetagents_langgraph_mcp_server --cov=regnetagents_langgraph_workflow --cov-report=term-missing`
6. Submit a pull request

### Code Guidelines

- Python 3.10+ compatible
- Include tests for new functionality
- Ensure existing tests pass before submitting

### Adding Domain Agents

See [docs/ADDING_NEW_CELL_TYPES.md](docs/ADDING_NEW_CELL_TYPES.md) for guidance on extending the system with new cell types or domain agents.

## Questions

Open an issue for questions about contributing.
