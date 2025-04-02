# dock-prep Documentation

Welcome to the dock-prep documentation. This tool helps prepare protein structures for molecular docking with AutoDock Vina.

## Table of Contents

- [Installation](installation/index.md)
  - [Prerequisites](installation/prerequisites.md)
  - [Installation Guide](installation/guide.md)
  - [Verifying Installation](installation/verification.md)

- [Usage Guide](usage/index.md)
  - [Basic Usage](usage/basic.md)
  - [Advanced Options](usage/advanced.md)
  - [Output Files](usage/output.md)
  - [Troubleshooting](usage/troubleshooting.md)

- [API Reference](api/index.md)
  - [Core Functions](api/core.md)
  - [Utilities](api/utilities.md)

- [Developer Guide](development/index.md)
  - [Setup](development/setup.md)
  - [Contributing](development/contributing.md)
  - [Testing](development/testing.md)

## Quick Start

```bash
# Install the package
git clone https://github.com/ingcoder/pdbqt-converter.git
cd dock-prep
pip install -e .

# Run the dependency checker
dock-prep-check

# Process a protein
dock-prep --pdb_id 2pgh --verbose
```

## Need Help?

If you encounter any issues or have questions:
- Check the [Troubleshooting](usage/troubleshooting.md) guide
- [Open an issue](https://github.com/ingcoder/pdbqt-converter/issues) on GitHub

## License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/ingcoder/pdbqt-converter/blob/main/LICENSE) file for details. 