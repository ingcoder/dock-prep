# Installation Guide

This section covers how to install dock-prep and its dependencies.

## Quick Install

For those who want to get started quickly:

```bash
# Create and activate a conda environment
conda create -n docking-pipeline python=3.10 -y
conda activate docking-pipeline

# Install dependencies
conda install -c conda-forge numpy pdbfixer openmm biopython openbabel pdb2pqr -y

# Install dock-prep
git clone https://github.com/ingcoder/pdbqt-converter.git
cd dock-prep
pip install -e .

# Install required tools
cd scripts
./install_mgltools.sh
./install_molprobity.sh
```

## Detailed Installation Instructions

For more detailed installation instructions, please refer to:

- [Prerequisites](prerequisites.md) - System requirements and dependencies
- [Installation Guide](guide.md) - Step-by-step installation process
- [Verifying Installation](verification.md) - How to check if the installation was successful

## Troubleshooting Installation

If you encounter issues during installation:

1. Check the [Troubleshooting](../usage/troubleshooting.md) guide
2. Run the dependency checker: `dock-prep-check`
3. [Open an issue](https://github.com/ingcoder/pdbqt-converter/issues) on GitHub if your problem persists 