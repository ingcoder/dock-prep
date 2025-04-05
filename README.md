# dock-prep: Automated Protein Preparation Tool for Molecular Docking with AutoDock Vina

[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Python 3.6+](https://img.shields.io/badge/python-3.6+-blue.svg)](https://www.python.org/downloads/)
[![GitHub issues](https://img.shields.io/github/issues/ingcoder/pdbqt-converter.svg)](https://github.com/ingcoder/pdbqt-converter/issues)


📋 [Features](#key-features) | 🚀 [Quick Start](#quick-start) | 📖 [Tutorial](#tutorial) | ⚙️ [Installation](#installation) | 🛠️ [Usage](#usage) | 📝 [Documentation](#documentation) | 🤝 [Contributing](#contributing)

![Protein preparation workflow](https://github.com/user-attachments/assets/4527e2ff-c436-4e94-ab1c-bfb86f0a4fed)
**Keywords**: protein structure preparation, molecular docking, AutoDock Vina, PDBQT conversion, computational drug discovery, PDB file processing

## **dock-prep** is a powerful, user-friendly tool that automates the preparation of protein structures for molecular docking simulations with AutoDock Vina, streamlining PDB to PDBQT conversion for computational drug discovery. Designed for researchers to quickly setup docking without the complex steps. 

## Why use dock-prep?
Molecular docking is a critical computational technique in drug discovery and structural biology that predicts how proteins interact with potential drug compounds (ligands). Before running any docking simulation with tools like AutoDock Vina, protein structures from the Protein Data Bank (PDB) require extensive preparation - a process that is technically challenging and time-consuming, especially for newcomers.

**dock-prep** eliminates these barriers by:

1. **Reducing technical complexity** - What normally requires multiple specialized tools and domain expertise becomes a simple, one-line command
2. **Saving significant time** - Preparation that typically takes hours of manual work is completed in minutes
3. **Minimizing user errors** - Automated protocols ensure consistent, reproducible results without human-introduced mistakes
4. **Providing a standardized workflow** - Creates a uniform process that ensures all structures meet quality standards for reliable docking results
5. **Maintaining full traceability** - Generates intermediate files at each step, allowing for quality control and methodological transparency

## Common Use Cases

**dock-prep** excels in numerous research and drug discovery scenarios:

### 💊 Structure-Based Drug Design
Rapidly prepare multiple protein targets for virtual screening campaigns of compound libraries against therapeutic targets. Screen thousands of compounds with confidence that your protein preparation won't introduce inconsistencies.

### 🧬 Protein-Protein Interaction Studies
Extract specific binding interfaces between protein chains using the distance-based selection feature, allowing focused analysis of protein-protein interaction sites.

### 🔬 Academic Research
Enable students and researchers without computational chemistry backgrounds to quickly incorporate molecular docking into their experimental workflows without needing to master complex preparation techniques.

### 🦠 Infectious Disease Research
Efficiently prepare viral protein structures (like COVID-19 targets) that often contain challenging features such as missing residues and non-standard amino acids for rapid therapeutic development.

### 🤖 High-Throughput Virtual Screening
Integrate into automated pipelines to process hundreds of protein structures consistently for large-scale computational screening efforts, maintaining quality across the entire dataset.

### 🧪 Experimental Structure Refinement
Clean and optimize newly solved experimental structures (X-ray, Cryo-EM) that may contain artifacts or preparation issues before conducting molecular simulations.

## Key Features

**dock-prep** solves common protein preparation challenges by:

- **Automating structure cleaning** - Removes unwanted water molecules, ions, and co-crystallized ligands
- **Repairing missing atoms and residues** - Completes protein structures with accurate modeling
- **Optimizing hydrogen placement** - Adds hydrogens at proper protonation states
- **Resolving steric clashes** - Fixes unfavorable conformations and atomic overlaps
- **Selecting binding sites** - Extracts specific protein chains by distance cutoffs or chain IDs
- **Assigning proper charges** - Adds essential charge parameters needed for simulation
- **Converting to PDBQT format** - Generates ready-to-dock files for AutoDock Vina

By handling these technical steps automatically, dock-prep allows researchers to focus on scientific questions rather than the complexities of file preprocessing.

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
