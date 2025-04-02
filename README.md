# dock-prep - Docking File Preparation Tool

[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Python 3.6+](https://img.shields.io/badge/python-3.6+-blue.svg)](https://www.python.org/downloads/)
[![GitHub issues](https://img.shields.io/github/issues/ingcoder/pdbqt-converter.svg)](https://github.com/ingcoder/pdbqt-converter/issues)

A streamlined tool for preparing protein structures for molecular docking with AutoDock Vina.

## Overview

This tool automates protein structure preparation for molecular docking by:
- Cleaning and completing PDB structures
- Extracting specific protein chains by cut-off distance or provided chain ID (useful to trim large proteins to binding site)
- Adding hydrogen atoms and optimizing side chains
- Converting optimized PDB files to PDBQT format for AutoDock Vina

## Tutorial

[Run the tutorial in Colab](https://colab.research.google.com/drive/1WDyGSLmT-XjFkU1L3d-mtd0GoD7p8EEy?usp=sharing)

## Installation

### 1. Set up Python Environment
```bash
conda create -n docking-pipeline python=3.10 -y
conda activate docking-pipeline
conda install -c conda-forge numpy pdbfixer openmm biopython openbabel pdb2pqr -y
```

### 2. Install PDBQT Converter
```bash
git clone https://github.com/ingcoder/pdbqt-converter.git
cd dock-prep
pip install -e .
```

### 3. Install Required Tools
```bash
# Install MGLTools
cd dock-prep/scripts
./install_mgltools.sh

# Install MolProbity
./install_molprobity.sh
```

### 4. Verify Installation
After installation, you can verify that all dependencies are properly installed:
```bash
# Run the dependency checker
dock-prep-check
```

This will check that:
- All required Python packages are installed
- You're running in a conda environment
- External tools (OpenBabel, PDB2PQR) are on your PATH
- Configuration-based tools (MGLTools, MolProbity) are properly configured

## Usage

Run the converter with a PDB ID or file:
```bash
# Process entire protein (default behavior)
dock-prep --pdb_id 2pgh --file_name path/to/protein.pdb --verbose

# Process specific chains
dock-prep --pdb_id 2pgh --file_name path/to/protein.pdb --target_chains A,B

# Extract chains by distance from a reference chain in angstrom
dock-prep --pdb_id 2pgh --file_name path/to/protein.pdb reference_chain A --distance 5.0
```

## Output Files

The tool produces the following output files in the current directory:
- `[pdb_id]_prepared.pdb`: Cleaned PDB with complete residues and optimized hydrogens
- `[pdb_id]_prepared.pdbqt`: Final PDBQT file ready for AutoDock Vina docking
- `[pdb_id]_validation.txt`: Structure validation report from MolProbity (optional)

## Troubleshooting

### Common MGLTools Issues

- **ImportError: No module named MolKit**: Ensure PYTHONPATH includes MGLToolsPckgs directory
- **No output file**: Check for error messages, verify input file exists, check write permissions

### Dependency Issues
If you encounter errors related to missing dependencies, run the dependency checker:
```bash
dock-prep-check
```
This will help identify which tools or packages need to be installed or properly configured.

## Dependencies

This tool relies on:
- **MGLTools**: For PDB to PDBQT conversion
- **MolProbity**: For structure validation and hydrogen placement
- **OpenBabel**: For file format conversion (`obabel`)
- **PDB2PQR**: For protein protonation (`pdb2pqr30`)

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. Check out our [Contributing Guidelines](CONTRIBUTING.md) for more details.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use this tool in your research, please cite:

```
Barbosa-Farias, I. (2023). dock-prep: A streamlined tool for preparing protein structures for molecular docking. 
GitHub repository: https://github.com/ingcoder/dock-prep
```

## Acknowledgments

- Thanks to all the developers of MGLTools, MolProbity, OpenBabel, and PDB2PQR
- Special thanks to contributors and users who have provided valuable feedback
