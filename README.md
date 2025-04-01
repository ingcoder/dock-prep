# PDBQT Converter

A streamlined tool for preparing protein structures for molecular docking with AutoDock Vina.

GoogleColab Tutorial: https://colab.research.google.com/drive/1WDyGSLmT-XjFkU1L3d-mtd0GoD7p8EEy?usp=sharing

## Overview

This tool automates protein structure preparation for molecular docking by:
- Cleaning and completing PDB structures
- Extracting specific protein chains by cut-off distance or provided chain ID (useful to trim large proteins to binding site)
- Adding hydrogen atoms and optimizing side chains
- Converting optimized PDB files to PDBQT format for AutoDock Vina

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
cd pdbqt-converter
pip install -e .
```

### 3. Install Required Tools
```bash
# Install MGLTools
cd pdbqt-converter/scripts
./install_mgltools.sh

# Install MolProbity
./install_molprobity.sh
```

## Usage

Run the converter with a PDB ID or file:
```bash
# Process entire protein (default behavior)
pdbqt-converter --pdb_id 2pgh --file_name path/to/protein.pdb --verbose

# Process specific chains
pdbqt-converter --pdb_id 2pgh --file_name path/to/protein.pdb --target_chains A,B

# Extract chains by distance from a reference chain in angstrom
pdbqt-converter --pdb_id 2pgh --file_name path/to/protein.pdb reference_chain A --distance 10.0
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

## Dependencies

This tool relies on:
- **MGLTools**: For PDB to PDBQT conversion
- **MolProbity**: For structure validation and hydrogen placement
- **OpenBabel**: For file format conversion (`obabel`)
- **PDB2PQR**: For protein protonation (`pdb2pqr30`)
