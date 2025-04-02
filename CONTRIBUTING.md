# Contributing to dock-prep

Thank you for considering contributing to dock-prep! This document provides guidelines and instructions to help you get started.

## Code of Conduct

By participating in this project, you agree to abide by our [Code of Conduct](CODE_OF_CONDUCT.md).

## How Can I Contribute?

### Reporting Bugs

Bug reports help us improve the tool. When reporting bugs:

1. **Check existing issues** to see if the problem has already been reported
2. **Use the bug report template** when creating a new issue
3. **Include detailed information** about your environment and steps to reproduce

### Suggesting Enhancements

Have ideas for improving dock-prep? We'd love to hear them!

1. **Check existing issues** to avoid duplicates
2. **Use the feature request template** when creating a new issue
3. **Be specific** about the enhancement and how it would benefit users

### Pull Requests

We welcome code contributions through pull requests:

1. **Fork the repository** and create your branch from `main`
2. **Install development dependencies** with `pip install -e ".[dev]"`
3. **Make your changes** following the coding standards
4. **Add tests** for any new features or bug fixes
5. **Run the test suite** to ensure your changes don't break existing functionality
6. **Update documentation** as needed
7. **Submit your pull request** with a clear description of the changes

## Development Setup

### Setting Up Your Environment

```bash
# Clone the repository
git clone https://github.com/YOUR_USERNAME/pdbqt-converter.git
cd dock-prep

# Create and activate a virtual environment
conda create -n docking-pipeline python=3.10 -y
conda activate docking-pipeline

# Install the package in development mode
pip install -e ".[dev]"

# Install required tools
cd scripts
./install_mgltools.sh
./install_molprobity.sh
```

### Running Tests

```bash
pytest
```

## Coding Standards

- Follow PEP 8 style guidelines
- Write docstrings for all functions, classes, and modules
- Comment complex code sections
- Keep functions focused and modular

## Commit Messages

- Use the present tense ("Add feature" not "Added feature")
- Use the imperative mood ("Move cursor to..." not "Moves cursor to...")
- Reference issues and pull requests where appropriate

## Pull Request Process

1. Update the README.md or documentation with details of changes if appropriate
2. Update the tests to reflect your changes
3. The PR should work for Python 3.6 and above
4. Your PR needs approval from at least one maintainer

Thank you for contributing to dock-prep! 