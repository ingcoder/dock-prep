# Using dock-prep

This section covers how to use the dock-prep tool for preparing protein structures for molecular docking.

## Basic Usage

The simplest way to use dock-prep is with a PDB ID:

```bash
dock-prep --pdb_id 2pgh --verbose
```

This will:
1. Download the PDB file
2. Clean the structure
3. Add hydrogens
4. Convert to PDBQT format

For more details, see [Basic Usage](basic.md).

## Command-Line Options

dock-prep provides various command-line options to customize the preparation process:

```bash
# Process specific chains
dock-prep --pdb_id 2pgh --target_chains A,B

# Extract chains by distance from a reference chain
dock-prep --pdb_id 2pgh --reference_chain A --distance 5.0

# Use a local PDB file
dock-prep --file_name path/to/protein.pdb
```

For a complete list of options, see [Advanced Options](advanced.md).

## Output Files

The tool produces several output files:
- `[pdb_id]_prepared.pdb`: Cleaned PDB with hydrogens
- `[pdb_id]_prepared.pdbqt`: Final PDBQT file for docking
- `[pdb_id]_validation.txt`: Validation report (if requested)

See [Output Files](output.md) for details on interpreting these files.

## Troubleshooting

If you encounter issues, check the [Troubleshooting](troubleshooting.md) guide. 