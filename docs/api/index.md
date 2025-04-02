# API Reference

This section documents the dock-prep API for developers who want to integrate the tool into their own Python scripts or extend its functionality.

## Module Structure

The dock-prep package is organized into the following modules:

- `dock_prep.run`: Entry point and command-line interface
- `dock_prep.pdb_prep`: PDB file preparation and cleaning
- `dock_prep.pdbqt_converter`: Conversion from PDB to PDBQT format
- `dock_prep.validation`: Structure validation utilities
- `dock_prep.utils`: Helper functions and utilities
- `dock_prep.check_dependencies`: Dependency checking functionality

## Programmatic Usage

You can use dock-prep programmatically in your Python scripts:

```python
from dock_prep.pdb_prep import clean_pdb
from dock_prep.pdbqt_converter import convert_to_pdbqt

# Clean a PDB file
cleaned_pdb_path = clean_pdb("input.pdb", add_hydrogens=True)

# Convert to PDBQT
pdbqt_path = convert_to_pdbqt(cleaned_pdb_path)

print(f"Created PDBQT file at: {pdbqt_path}")
```

## Core Functions

See [Core Functions](core.md) for documentation on the main functions for protein preparation.

## Utility Functions

See [Utilities](utilities.md) for documentation on helper functions and utilities.

## Extension Points

dock-prep is designed to be extensible. Common extension points include:

1. Adding new structure validation methods
2. Implementing custom PDB cleaning protocols
3. Adding support for new formats or docking software

Refer to the [Developer Guide](../development/index.md) for more information on extending the package. 