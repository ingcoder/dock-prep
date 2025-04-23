# dock-prep: Automated Protein Preparation Tool for Molecular Docking with AutoDock Vina

[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Python 3.6+](https://img.shields.io/badge/python-3.6+-blue.svg)](https://www.python.org/downloads/)
[![GitHub issues](https://img.shields.io/github/issues/ingcoder/pdbqt-converter.svg)](https://github.com/ingcoder/pdbqt-converter/issues)


📋 [Features](#key-features) | 🚀 [Quick Start](#quick-start) | 📖 [Tutorial](#tutorial) | ⚙️ [Installation](#installation) | 🛠️ [Usage](#usage) | 📝 [Documentation](#documentation) | 🤝 [Contributing](#contributing)

![Protein preparation workflow](https://github.com/user-attachments/assets/4527e2ff-c436-4e94-ab1c-bfb86f0a4fed)
**Keywords**: protein structure preparation, molecular docking, AutoDock Vina, PDBQT conversion, computational drug discovery, PDB file processing

## What is dock-prep?
**dock-prep** is a powerful, user-friendly tool that automates the preparation of protein structures for molecular docking with AutoDock Vina, streamlining PDB to PDBQT creation for computational drug discovery. Designed for researchers to convert PDBs as they come from the Protein Databank to ready to use PDBQT files for AutoDock Vina docking, in a single line of code. 

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
Process protein structures consistently for large-scale screening with standardized protocols. 

## Quick Start

```bash
# Create and activate conda environment
conda create -n docking python=3.10 -y && activate docking

# Install dependencies and dock-prep
conda install -c conda-forge numpy pdbfixer openmm biopython openbabel pdb2pqr -y
git clone https://github.com/ingcoder/dock-prep.git
pip install -e dock-prep

# Install external tools
chmod +x dock-prep/scripts/*.sh
./dock-prep/scripts/install_mgltools.sh
./dock-prep/scripts/install_molprobity.sh #optional, but recommended

# Prepare a protein from PDB ID
dock-prep --input_file dock-prep/dock_prep/examples/1n6d.pdb --reference_atom_chains H --cutoff 2.0 --verbose
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

## Usage (Linux or MacOS)

1️⃣ **Download pdb** from the Protein Data Bank.

2️⃣ **Move pdb** into your project folder.

3️⃣ **Activate conda** environment (if not already active):   
```bash 
conda activate dock-prep
```

4️⃣ **Run dock-prep** with one of the dock-prep commands shown below, e.g. 
```bash 
dock-prep --file_input MyProjectFolder/1abc.pdb --verbose
```

5️⃣ **Check results** The processed file are in the automatically created **results/** folder inside your project directory.

Your file structure will look like this. 
```text
MyProjectFolder/
├── dock-prep/                      # Dock-Prep repo or package
├── 1abc.pdb                        # raw input structure
└── results/
    └── 1abc_structure_docking.pdbqt
 ```

## Usage (Colab Notebook)
If you want to run the colab notebook using your own pdb file: 
[Run the interactive tutorial in Google Colab](https://colab.research.google.com/drive/1WDyGSLmT-XjFkU1L3d-mtd0GoD7p8EEy?usp=sharing)

1️⃣ **Copy Notebbook** Open link above and copy notebooke with: File -> Save a copy in drive 

2️⃣ **Install dock-prep** Run all cells in installation section to install dock-prep and dependencies

3️⃣ **Download pdb** from the Protein Data Bank.

4️⃣ **Upload pdb** file to colab. Click the folder icon in the sidebar, then the ⬆️ upload button. The project folder in colab is called content/.
<img width="200" alt="Image" src="https://github.com/user-attachments/assets/4dc662f8-f628-4538-b801-842456a3bdfa" />

5️⃣ **Run dock-prep** replace the name of the pdb file with your filename and run the cell containing dock-prep --file_input content/1abc.pdb --verbose --skip_molprobity. 
Replace filename with your filename.

6️⃣ **Check results** in the automatically created **results/** folder inside the content/ directory.

Your file structure will look like this. 
```text
content/
├── dock-prep/                      # Dock-Prep repo or package
├── 1abc.pdb                        # raw input structure
└── results/
    └── 1abc_structure_docking.pdbqt
```


### Basic Commands

Run the converter with a PDB ID or file:
```bash
# Process entire protein (default behavior, works for small proteins)
dock-prep --file_input path/to/1abc.pdb --verbose

# Process specific chains
dock-prep --file_input path/to/1abc.pdb --include_chains A,B --verbose

# Extract chains by distance from a reference peptide chain in angstrom (5 Angstrom by default)
dock-prep --file_input path/to/1abc.pdb reference_atom_chains H --cutoff 2.0 --verbose

# Extract chains by distance from a reference small molecule hetatom chain in angstrom (5 Angstrom by default)
dock-prep --file_input path/to/1abc.pdb reference_hetatm_chains H --cutoff 2.0 --verbose
```

### Output Files 

The tool generates a series of progressively refined files that document each step in the protein preparation pipeline:

| File | Description | Purpose in Workflow |
|------|-------------|---------------------|
| 📄`_structure_cleaned.pdb` | Initial cleaned structure | Removes HETATM records (waters, ligands, ions) and prepares the protein for structural completion |
| 📄`_structure_completed_final.pdb` | Structure with modeled residues | Fills in missing atoms and residues to create a complete protein model |
| 📄`_structure_flipped_h_final.pdb` | Optimized with hydrogens | Contains MolProbity-optimized hydrogen positions and corrected side-chain orientations |
| 📄`_structure_protonated.pqr` | Protonated structure | Includes atomic radii and charge parameters from PDB2PQR required for electrostatics |
| 📄`_structure_docking.pdbqt` | Final docking-ready file | **Primary output file** with all parameters needed for AutoDock Vina docking simulations |

> **Note**: The `docking.pdbqt` file is the primary output that should be used for docking simulations with AutoDock Vina. Intermediate files are preserved to allow inspection of each preparation step.

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
