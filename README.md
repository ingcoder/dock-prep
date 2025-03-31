# Protein Docking Pipeline

A comprehensive suite of scripts and tools for protein structure preparation and molecular docking.

## Table of Contents
- [Overview](#overview)
- [Installation](#installation)
  - [Docking Pipeline Environment](#docking-pipeline-environment)
  - [MGLTools Environment](#mgltools-environment)
  - [MolProbity Tools](#setting-up-molprobity-tools)
- [Usage](#usage)
  - [Basic Workflow](#basic-workflow)
  - [File Preparation](#file-preparation)
  - [PDB to PDBQT Conversion](#pdb-to-pdbqt-conversion)
- [Tool Documentation](#tool-documentation)
  - [MolProbity Tools](#molprobity-tools)
  - [MGLTools](#mgltools)
- [Troubleshooting](#troubleshooting)
- [Notes](#notes)

## Overview

This repository contains scripts and workflows for protein structure preparation and molecular docking, integrating several established tools:

- **Structure Preparation**: Clean PDBs, complete missing residues, restore chain IDs, and add protons
- **Structure Validation**: Assess and improve structural quality with MolProbity
- **Docking Preparation**: Convert files to PDBQT format for docking software

## Installation

### Docking Pipeline Environment

The main pipeline runs in a Python 3.10 environment with several specialized bioinformatics packages.

```bash
# Create and activate the conda environment
conda create -n docking-pipeline python=3.10 -y
conda activate docking-pipeline

# Install required packages
conda install -c conda-forge numpy pdbfixer openmm biopython openbabel pdb2pqr jupyter -y
```

### Clone Github repository 
``bash
git clone https://github.com/ingcoder/pdbqt-converter.git
cd pdbqt-converter
```

### Setting Up Dependencies
#1 MolProbity Tools
MolProbity provides tools for structure validation and hydrogen placement.


```bash
# Clone and install MolProbity
git clone https://github.com/rlabduke/MolProbity.git
cd MolProbity

# Make the installation script executable
chmod +x ./install_via_bootstrap.sh 4

# Run the setup script
./install_via_bootstrap.sh 4

# Apply the changes
source ~/.zshrc

# Add MolProbity bin to PATH
export PATH=\$PATH:/Users/ingrid/PycharmProjects/jupyternotebook/MolProbity/molprobity/bin >> ~/.zshrc
```


### MGLTools Environment

The PDB to PDBQT conversion requires MGLTools, which needs Python 2.7.

```bash
# Create Python 2.7 environment for MGLTools
conda create -n mgltools python=2.7 -y
conda activate mgltools
```

# Installation
https://ccsb.scripps.edu/mgltools/downloads/
Instructions: https://yzhang.hpc.nyu.edu/DeltaVina/tutorial.html

Download to working directory. mgltools_1.5.7_MacOS-X.tar.gz (tarball installer 85Mb)

```bash
tar -xvzf mgltools_x86_64Linux2_1.5.6.tar.gz
cd mgltools_x86_64Linux2_1.5.6/
./install.sh
```

After installing MGLTools, set up the environment variables:

```bash
# Set MGLTools paths (modify with your actual installation path)
export MGLTOOLS_HOME=/path/to/mgltools_1.5.7_MacOS-X
export PYTHONPATH=$MGLTOOLS_HOME/MGLToolsPckgs:$PYTHONPATH
export PATH=$MGLTOOLS_HOME/bin:$PATH
```

Example with actual path:
```bash
export PYTHONPATH=/Users/ingrid/PycharmProjects/jupyternotebook/mgltools_1.5.7_MacOS-X/MGLToolsPckgs:$PYTHONPATH
```



## Usage

### Basic Workflow

```bash
# Activate the environment
conda activate docking-pipeline  # or use the 'dp' alias

# Run the main preparation script
python file_handler.py
```

### File Preparation

The `file_handler.py` script performs these operations:
- Structure cleaning and binding site extraction
- Missing residue completion
- Chain ID restoration
- Structure optimization with MolProbity
- Protonation using PDB2PQR
- File format conversion (PQR to PDB)

### PDB to PDBQT Conversion

For molecular docking, PDB files need to be converted to PDBQT format using MGLTools.

#### Converting Receptor Files

```bash
# Activate MGLTools environment
conda activate mgltools

# Convert protein receptor to PDBQT
$MGLTOOLS_HOME/bin/pythonsh $MGLTOOLS_HOME/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py \
  -r protein.pdb \
  -o protein.pdbqt \
  -A hydrogens
```

Example with absolute paths:
```bash
/Users/ingrid/PycharmProjects/jupyternotebook/mgltools_1.5.7_MacOS-X/bin/pythonsh \
/Users/ingrid/PycharmProjects/jupyternotebook/mgltools_1.5.7_MacOS-X/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py \
-r your_file.pdb -o output.pdbqt -A hydrogens
```

#### Converting Ligand Files

```bash
$MGLTOOLS_HOME/bin/pythonsh $MGLTOOLS_HOME/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py \
  -l ligand.pdb \
  -o ligand.pdbqt \
  -A hydrogens
```

#### Common Options for PDB to PDBQT Conversion
- `-r`/`-l`: Input PDB file (receptor/ligand)
- `-o`: Output PDBQT file
- `-A`: Add specific atoms (hydrogens)
- `-C`: Add charges
- `-p`: Preserve input charges

## Tool Documentation

### MolProbity Tools

MolProbity includes utilities for structural validation and hydrogen placement:

#### Reduce: Hydrogen Addition and Optimization

```bash


# Add hydrogens without flips
reduce -NOFLIP input.pdb > output_with_H.pdb 

# Remove hydrogens
reduce -TRIM input.pdb > output_no_H.pdb
```

#### Probe: Steric Clash Detection

```bash
# Generate all-atom contacts
probe -u -self input.pdb > contacts.probe

# Generate clash score
probe -4 -u -self input.pdb > clashes.probe
```

### MGLTools

MGLTools provides utilities for preparing files for molecular docking:

#### Example 1: Convert Protein
```bash
/path/to/mgltools/bin/pythonsh \
/path/to/mgltools/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py \
-r 1hsg.pdb -o 1hsg.pdbqt -A hydrogens
```

#### Example 2: Convert Ligand
```bash
/path/to/mgltools/bin/pythonsh \
/path/to/mgltools/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py \
-l indinavir.pdb -o indinavir.pdbqt -A hydrogens
```

## Troubleshooting

### MGLTools Issues

#### ImportError: No module named MolKit
- Ensure PYTHONPATH includes the MGLToolsPckgs directory
- Use pythonsh from MGLTools, not system Python
- Verify paths are correct for your installation

#### No output file generated
- Check for error messages in terminal
- Verify input PDB file exists and is valid
- Check write permissions in output directory

### OpenBabel Issues
- The OpenBabel conversion sometimes times out, in which case the script falls back to an internal Python implementation

## Notes

- The MolProbity tools are installed in the MolProbity directory
- This pipeline requires external command-line tools that are included in the conda environment setup:
  - `obabel` (from OpenBabel package): For file format conversion
  - `pdb2pqr30` (from PDB2PQR package): For protein protonation