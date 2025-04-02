# Developer Guide

This section provides information for developers who want to contribute to dock-prep or extend its functionality.

## Development Setup

To set up a development environment:

```bash
# Clone the repository
git clone https://github.com/YOUR_USERNAME/pdbqt-converter.git
cd dock-prep

# Create and activate a conda environment
conda create -n docking-pipeline-dev python=3.10 -y
conda activate docking-pipeline-dev

# Install in development mode with development dependencies
pip install -e ".[dev]"

# Install required external tools
cd scripts
./install_mgltools.sh
./install_molprobity.sh
```

For more details, see [Setup](setup.md).

## Running Tests

We use pytest for testing:

```bash
# Run all tests
pytest

# Run tests with coverage
pytest --cov=dock_prep tests/

# Run a specific test file
pytest tests/test_pdb_prep.py
```

For more information on testing, see [Testing](testing.md).

## Code Style

We follow PEP 8 style guidelines. You can check and format your code using:

```bash
# Check code style
flake8 dock_prep/

# Format code
black dock_prep/
```


## Contributing

Please read our [Contributing Guidelines](../../CONTRIBUTING.md) for details on the process for submitting pull requests.

## Release Process

For maintainers, the release process is:

1. Update version in `setup.py`
2. Update CHANGELOG.md
3. Create a new GitHub release
4. Publish to PyPI (if applicable)

See [Release Process](release.md) for more details. 