# dock-prep: Automated Protein Preparation Tool for Molecular Docking with AutoDock Vina

[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Python 3.6+](https://img.shields.io/badge/python-3.6+-blue.svg)](https://www.python.org/downloads/)
[![GitHub issues](https://img.shields.io/github/issues/ingcoder/pdbqt-converter.svg)](https://github.com/ingcoder/pdbqt-converter/issues)


📋 [Features](#key-features) | 🚀 [Quick Start](#quick-start) | 📖 [Tutorial](#tutorial) | ⚙️ [Installation](#installation) | 🛠️ [Usage](#usage) | 📝 [Documentation](#documentation) | 🤝 [Contributing](#contributing)

![Protein preparation workflow](https://github.com/user-attachments/assets/4527e2ff-c436-4e94-ab1c-bfb86f0a4fed)
**Keywords**: protein structure preparation, molecular docking, AutoDock Vina, PDBQT conversion, computational drug discovery, PDB file processing

## What is dock-prep?
**dock-prep** is a powerful, user-friendly tool that automates the preparation of protein structures for molecular docking with AutoDock Vina, streamlining PDB to PDBQT creation for computational drug discovery. Designed for researchers to quickly setup docking without the complex steps. 

## Why use dock-prep?
**dock-prep** handles the entire pipeline of file preperation, reducing manual errors, ensuring consistency and saving time. 

#### Key Features & Benefits

| Feature                      | What It Does                                         | Why It Matters                                                               |
|---------------------------  |------------------------------------------------------|------------------------------------------------------------------------------|
| ✅ **Structure Cleaning**   | Removes waters, ions, and ligands                    | Avoids docking to irrelevant or non-biological parts  |
| ✅ **Gap filling**          |   Completes missing atoms and residues                | Docking tools require complete structures   |
| ✅ **Hydrogen Addition**    | Adds hydrogens with protonation at pH             | Ensures accurate hydrogen bonding prediction                 |
| ✅ **Clash Resolution**    | Fixes unfavorable sidechain conformations            | Reduces steric clashes that could disrupt key interactions              |
| ✅ **Site Selection**      | Extracts chains by chain IDs or distance to ligand                 | Focuses on biologically meaningful interaction regions            |
| ✅ **Charge Assignment**        | Assigns atomic charges and radii                   | Enables MD simulations requiring charge information   |
| ✅ **PDBQT File Conversion**     | Generates PDBQT files                  | Provides required PDBQT format for AutoDock Vina             |


## Common Use Cases

**dock-prep** excels in numerous research and drug discovery scenarios:

#### 💊 Structure-Based Drug Design
Process experimental structures into docking-ready models with optimized parameters for accurate virtual screening across diverse protein families.
#### 🧬 Protein Trimming for Focused Docking
Extract only relevant binding pockets through manual chain selection or distance-based trimming to improve docking accuracy and computational efficiency.
#### 🔬🏢 Research Across Academia and Industry
Simplify molecular docking for academic teaching and biotech R&D teams while reducing computational costs and accelerating drug discovery timelines.
#### 🤖 High-Throughput Virtual Screening Pipelines
Process hundreds of protein structures consistently for large-scale screening with standardized protocols that integrate with existing high-performance computing environments.

## Quick Start

```bash
# Install dependencies with conda
conda create -n docking python=3.10 -y && conda activate docking
conda install -c conda-forge numpy pdbfixer openmm biopython openbabel pdb2pqr -y

# Install dock-prep
git clone https://github.com/ingcoder/dock-prep.git
pip install -e dock-prep 

# Install external tools
chmod +x dock-prep/scripts/*.sh
./dock-prep/scripts/install_mgltools.sh
./dock-prep/scripts/install_molprobity.sh

# Prepare a protein from PDB ID
dock-prep --pdb_id 2pgh --input_file dock-prep/dock_prep/examples/2pgh_original.pdb --verbose
```

## Tutorial

Follow our comprehensive tutorial to learn how dock-prep can integrate into your molecular docking workflow:

[Run the interactive tutorial in Google Colab](https://colab.research.google.com/drive/1WDyGSLmT-XjFkU1L3d-mtd0GoD7p8EEy?usp=sharing)

## Installation

### 1. Set up Python Environment
```bash
conda create -n docking-pipeline python=3.10 -y
conda activate docking-pipeline
conda install -c conda-forge numpy pdbfixer openmm biopython openbabel pdb2pqr -y
```

### 2. Install dock-prep 
```bash
git clone https://github.com/ingcoder/dock-prep.git
pip install -e dock-prep 
```

### 3. Install Required Tools
```bash
# Install MGLTools
cd dock-prep/scripts
chmod +x install_mgltools.sh  # Ensure script has executable permissions
./install_mgltools.sh

# Install MolProbity
chmod +x install_molprobity.sh  # Ensure script has executable permissions
./install_molprobity.sh
```

> **Important Note**: If you encounter "permission denied" errors when running the scripts, you need to manually set executable permissions using the `chmod +x script_name.sh` command. The scripts include self-fixing permission code, but this only works if the script can be executed in the first place.

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

### Basic Commands

Run the converter with a PDB ID or file:
```bash
# Process entire protein (default behavior)
dock-prep --pdb_id PDBID --file_name path/to/protein.pdb --verbose

# Process specific chains
dock-prep --pdb_id PDBID --file_name path/to/protein.pdb --target_chains A,B

# Extract chains by distance from a reference chain in angstrom
dock-prep --pdb_id PDBID --file_name path/to/protein.pdb reference_chain A --distance 5.0
```

### Output Files 

The tool generates a series of progressively refined files that document each step in the protein preparation pipeline:

| File | Description | Purpose in Workflow |
|------|-------------|---------------------|
| `[pdb_id]_structure_cleaned.pdb` | Initial cleaned structure | Removes HETATM records (waters, ligands, ions) and prepares the protein for structural completion |
| `[pdb_id]_structure_completed_final.pdb` | Structure with modeled residues | Fills in missing atoms and residues to create a complete protein model |
| `[pdb_id]_structure_flipped_h_final.pdb` | Optimized with hydrogens | Contains MolProbity-optimized hydrogen positions and corrected side-chain orientations |
| `[pdb_id]_structure_protonated.pqr` | Protonated structure | Includes atomic radii and charge parameters from PDB2PQR required for electrostatics |
| `[pdb_id]_structure_docking.pdbqt` | Final docking-ready file | **Primary output file** with all parameters needed for AutoDock Vina docking simulations |

> **Note**: The `[pdb_id]_refined_docking.pdbqt` file is the primary output that should be used for docking simulations with AutoDock Vina. Intermediate files are preserved to allow inspection of each preparation step.

#### File Formats Explained

- **PDB**: Standard Protein Data Bank format containing atomic coordinates
- **PQR**: Modified PDB format that includes charge (Q) and radius (R) parameters
- **PDBQT**: Extended PDB format with partial charges (Q) and atom types (T) required by AutoDock Vina

For advanced users who want to customize the preparation process, these intermediate files can be modified before continuing to the next processing step using the `--input_file` parameter.

## Documentation

### Troubleshooting

#### Common MGLTools Issues

- **ImportError: No module named MolKit**: Ensure PYTHONPATH includes MGLToolsPckgs directory
- **No output file**: Check for error messages, verify input file exists, check write permissions

#### Dependency Issues
If you encounter errors related to missing dependencies, run the dependency checker:
```bash
dock-prep-check
```
This will help identify which tools or packages need to be installed or properly configured.

### Dependencies

This tool relies on:
- **MGLTools**: For PDB to PDBQT conversion
- **MolProbity (optional)**: For structure validation and hydrogen placement
- **OpenBabel**: For file format conversion (`obabel`)
- **PDB2PQR**: For protein protonation (`pdb2pqr30`)

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. Check out our [Contributing Guidelines](CONTRIBUTING.md) for more details.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use this tool in your research, please cite:

```
Barbosa-Farias, I. (2025). dock-prep: A streamlined tool for preparing protein structures for molecular docking. 
GitHub repository: https://github.com/ingcoder/dock-prep
```

## Acknowledgments

- Thanks to all the developers of MGLTools, MolProbity, OpenBabel, and PDB2PQR
- Special thanks to contributors and users who have provided valuable feedback
