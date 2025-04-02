# Contributing to dock-prep

Thank you for considering contributing to dock-prep!

## Code of Conduct

Contributors must follow our [Code of Conduct](CODE_OF_CONDUCT.md).

## How to Contribute

### Bugs & Features
- Check existing issues first
- Use templates for bug reports and feature requests
- Provide specific details

### Pull Requests
1. Fork from `main` branch
2. Install dev dependencies: `pip install -e ".[dev]"`
3. Follow coding standards & add tests
4. Update documentation as needed

## Development

```bash
# Setup
git clone https://github.com/YOUR_USERNAME/pdbqt-converter.git
cd dock-prep
conda create -n docking-pipeline python=3.10 -y && conda activate docking-pipeline
pip install -e ".[dev]"

# Install tools
cd scripts && ./install_mgltools.sh && ./install_molprobity.sh

# Test
pytest
```

## Standards
- Follow PEP 8
- Write docstrings
- Use present tense & imperative mood in commits

## PR Requirements
- Update docs & tests
- Support Python 3.6+
- Get approval from a maintainer

Thank you for contributing! 