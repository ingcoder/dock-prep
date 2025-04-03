"""
Structure Handler for Molecular Docking Preparation

This script prepares protein structures for molecular docking by performing a comprehensive 
workflow including structure cleaning, binding site extraction, missing structure completion,
and format conversion for docking software.

The workflow consists of these main steps:
1. Loading and cleaning the original PDB structure
2. Extracting the binding site based on specified criteria
3. Identifying and completing missing structure elements
4. Optimizing hydrogen positions for accurate docking 
5. Converting to PDBQT format for docking

Dependencies:
    - numpy, pdbfixer, openmm, biopython
    - External tools: MolProbity, MGLTools, OpenBabel

Author: Ingrid Barbosa-Farias
Date: 2025-03-27
"""

import os
import sys
from dock_prep.subprocess_handler import run_program

try:
    import numpy as np
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
    from openmm import unit
    from Bio import PDB
except ImportError as e:
    print(f"ERROR: Missing dependency - {e}")
    print("\nPlease install the required dependencies:")
    print("  - numpy: pip install numpy")
    print("  - pdbfixer: conda install -c conda-forge pdbfixer")
    print("  - openmm: conda install -c conda-forge openmm")
    print("  - biopython: pip install biopython")
    print("\nPlease install the missing dependencies and try again.")
    sys.exit(1)

#==============================================================================
#                            UTILITY FUNCTIONS
#==============================================================================

def add_separator(message):
    """Prints a separator with a centered message."""
    print("\n" + "="*80)
    print(message)
    print("="*80)

#==============================================================================
#                            FILE OPERATIONS
#==============================================================================

def extract_chains_to_pdb(input_filepath, output_filepath, target_chains=None):
    """Extracts specified chains from a PDB file and saves to a new file."""
    with open(input_filepath, 'r') as fin, open(output_filepath, 'w') as fout:
        for line in fin:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                chain_id = line[21]
                if chain_id in target_chains:
                    fout.write(line)
            elif line.startswith("TER") and len(line) > 21:
                chain_id = line[21]
                if chain_id in target_chains:
                    fout.write(line)
            else:
                fout.write(line)
    
    rel_path = os.path.relpath(output_filepath)
    print(f"✅ Filtered PDB saved to {rel_path} containing chains {target_chains}")

def save_structure_to_pdb(pdb_structure, output_filepath, verbose=True):
    """Saves PDBFixer structure to a PDB file."""
    with open(output_filepath, 'w') as f:
        PDBFile.writeFile(pdb_structure.topology, pdb_structure.positions, f)
    if verbose:
        rel_path = os.path.relpath(output_filepath)
        print(f"\n✅ Complete structure saved to {rel_path}")

def restore_original_chain_ids(input_filepath, output_filepath, target_chains, verbose=True):
    """
    Restores original chain IDs in a PDB file created by PDBFixer.
    
    PDBFixer sometimes changes the chain IDs during processing. This function
    maps the new chain IDs back to the original ones specified by the user.
    

    """
    # Parse the input PDB file using Biopython's PDB parser
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("orig", input_filepath)
    
    # Get all chain IDs from the input file by iterating through all models and chains
    input_file_chain_ids = [chain.id for model in structure for chain in model]
    
    # If no target chains were specified, use all chains from the input file
    if target_chains is None:
        target_chains = input_file_chain_ids
    
    # Create a mapping between the chains in the input file and the target chains
    # This assumes chains are in the same order in both lists (PDBFixer typically preserves order)
    chain_mapping = {orig_chain: new_chain for orig_chain, new_chain in zip(input_file_chain_ids, target_chains)}
    verbose and print(f"• Chain mapping: {chain_mapping}")
    
    with open(input_filepath, 'r') as fin, open(output_filepath, 'w') as fout:
         for line in fin:
             # Process only ATOM and HETATM lines which contain atom coordinates
             if line.startswith("ATOM") or line.startswith("HETATM"):
                 # In PDB format, the chain ID is at position 21 (0-indexed)
                 new_chain = line[21]
                 # Only modify the chain ID if it has a mapping
                 if new_chain in chain_mapping:
                     # Get the target chain ID from the mapping
                     fixed_chain = chain_mapping[new_chain]
                     # Construct a new line by replacing just the chain ID character
                     # This preserves all other information in the line
                     line = line[:21] + fixed_chain + line[22:]
             # Write the line (either modified or original) to the output file
             fout.write(line)
    if verbose:
        rel_path = os.path.relpath(output_filepath)
        print(f"✅ Fixed PDB saved to {rel_path}")

#==============================================================================
#                     STRUCTURE LOADING AND CLEANING
#==============================================================================

def load_clean_structure(input_file, output_file, no_hetatm=True, verbose=True):
    """Loads PDB file, removes heteroatoms, and returns a PDBFixer object."""
    # Validate input file
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"File not found. {input_file}")
    if verbose:
        # Option 1: Get path relative to current working directory
        rel_path = os.path.relpath(input_file)
        print("\nPath to PDB file:", rel_path)
        
        # Other options (commented out):
        # Option 2: Just the filename without the directory path
        # print("• Filename:", os.path.basename(pdb_filepath))
        
        # Option 3: Relative to a specific directory
        # print("• Path relative to results folder:", os.path.relpath(pdb_filepath, 'results'))
        
        # Option 4: Shorten the path by showing ~ for home directory
        # home = os.path.expanduser("~")
        # if pdb_filepath.startswith(home):
        #     shortened_path = pdb_filepath.replace(home, "~", 1)
        #     print("• Path to PDB file:", shortened_path)
        # else:
        #     print("• Path to PDB file:", pdb_filepath)
    
    # # Prepare output filename
    # base, ext = os.path.splitext(pdb_filepath)
    # basename = os.path.basename(base)
    # print(f"• Base: {basename}")
    # print(os.getcwd)
    # print(f"• Ext: {ext}")
    # cleaned_filepath = f'{base}_clean{ext}'
    
    try:
        # Step 1: Read and filter PDB content
        with open(input_file, 'r') as f:
            lines = f.readlines()

        # Step 2: Process lines and filter HETATM records
        cleaned_lines = []
        for line in lines:
            if line.startswith("HETATM"):
                if no_hetatm:  
                    # Skip all HETATM lines
                    continue
                else:
                    # Skip water molecules only
                    water_molecules = ["HOH", "WAT", "H2O", "SOL", "DOD"]
                    if any(water in line for water in water_molecules):
                        continue
            cleaned_lines.append(line)
    
        # Step 3: Write cleaned content to file
        cleaned_pdb_text = "".join(cleaned_lines)
        with open(output_file, 'w') as f:
            f.write(cleaned_pdb_text)

        return load_structure_as_pdbfixer(input_file)
    except Exception as e:
        print(f"❌ Error processing PDB file: {e}")
        raise
        

def load_structure_as_pdbfixer(input_file):
        # Step 4: Create and return PDBFixer object
    try:
        fixer = PDBFixer(filename=input_file)
        return fixer
    except Exception as e:
        print(f"❌ Error processing PDB file: {e}")
        raise

def count_residues(fixer, verbose=True):
    """Counts residues and atoms in a PDBFixer object."""
    residue_count = 0
    atom_count = 0
    for chain in fixer.topology.chains():
        for residue in chain.residues():
            residue_count += 1
            for atom in residue.atoms():
                atom_count += 1
    
    if verbose:
        print("• Residue count:", residue_count)
        print("• Atom count:", atom_count)
    
    return residue_count, atom_count

#==============================================================================
#                        STRUCTURE REFINEMENT
#==============================================================================

def complete_missing_structure(pdb_fixer, missing_residues_dict=None, verbose=True):
    """Completes missing residues and atoms in a PDBFixer object."""
    if verbose:
        print("\n> Before structure completion:")
    
    residue_count_before, atom_count_before = count_residues(pdb_fixer, verbose=verbose)
    
    # Step 1: Identify missing residues
    pdb_fixer.findMissingResidues()
    
    # Step 2: Handle missing residues (automatic or from dictionary)
    if missing_residues_dict is None:
        auto_missing_residues = pdb_fixer.missingResidues
        if verbose:
            print(f"\n> Automatically detected missing residues: {len(auto_missing_residues)}")
            for chain in auto_missing_residues:
                print(f"• Chain {chain}: {len(auto_missing_residues[chain])} missing residue segments")
    else:
        # Use provided missing residues dictionary
        pdb_fixer.missingResidues = missing_residues_dict
        if verbose:
            print(f"\n> Sanity check: number of chains with missing residues: {len(pdb_fixer.missingResidues)} == {len(missing_residues_dict)}")
    
    # Step 3: Identify and report missing atoms and terminals
    pdb_fixer.findMissingAtoms()
    missing_atoms = pdb_fixer.missingAtoms
    missing_terminals = pdb_fixer.missingTerminals
    if verbose:
        print(f"> Number of residues with missing atoms found: {len(missing_atoms)} residues with missing atoms")
        print(f"> Number of chains with missing terminals found: {len(missing_terminals)} chains with missing terminals")
    
    # Step 4: Add missing atoms to the structure
    pdb_fixer.addMissingAtoms()
    # Hydrogens can be added here if needed:
    # pdb_fixer.addMissingHydrogens(pH=6.2)

    if verbose:
        print("\n> After structure completion:")
    
    residue_count_after, atom_count_after = count_residues(pdb_fixer, verbose=verbose)

    return pdb_fixer, residue_count_before, atom_count_before, residue_count_after, atom_count_after

def get_missing_residues_by_chain(pdb_fixer, target_chains, verbose=True):
    """Gets missing residues in specified chains."""
    summaries = []
    missing_residues_in_target_chains = {}
    
    # Step 1: Create a mapping from chain indices (int) to chain objects
    chain_index_map = {}
    for chain in pdb_fixer.topology.chains():
        if target_chains is not None:
            if chain.id in target_chains:
                # verbose and print('target chains', target_chains, chain.id)
                # verbose and print(f"• Chain {chain.id} is in target chains")
                chain_index_map[chain.index] = chain
        else: # If no target chains are specified, include all chains
            chain_index_map[chain.index] = chain
    
    # Step 2: Collect missing residues for (target) chains
    for (chain_index, residue_index), residues in pdb_fixer.missingResidues.items():
        if chain_index in chain_index_map:
            # Skip residue index 0 (often used for special cases)
            if residue_index != 0:  
                chain = chain_index_map[chain_index]
                summary = {
                    'chain_id': chain.id,
                    'chain_index': chain.index,
                    'residue_index': residue_index,
                    'missing_residues': residues,
                    'missing_residues_count': len(residues)
                }
                summaries.append(summary)
                missing_residues_in_target_chains[(chain_index, residue_index)] = residues

    # Step 3: Print summary report of missing residues
    if verbose:
        print(f"\nMissing residues by target chains:")
        print(f"\n{'Chain ID':<10} {'Missing residues':<15} {'Position':<10}")
        print(f"{'-'*10} {'-'*15} {'-'*10}")
        for summary in summaries:
            print(f"{summary['chain_id']:<10} {summary['missing_residues_count']:<15} {summary['residue_index']:<10}")

    return missing_residues_in_target_chains

#==============================================================================
#                      CHAIN AND RESIDUE ANALYSIS
#==============================================================================

def get_selected_chains(input_file, selection, fixer_object, cutoff):
    selected_chains = None

        # Map the original chain IDs to the potentially new IDs after cleaning
    chain_ids = map_orig_to_new_chain(input_file, selection)
    # Identify atoms belonging to the specified HETATM chains
    ligand_atoms = get_ligand_atoms(fixer_object, chain_ids) 
    # Find protein residues within the cutoff distance of the HETATM atoms
    residues_near_ligand = get_residues_near_ligand(
        fixer_object, ligand_atoms, distance_threshold=cutoff
    )
    # Determine which chains contain the identified nearby residues
    selected_chains = get_chains_near_ligand(fixer_object, residues_near_ligand, chain_ids)

    return selected_chains


def map_orig_to_new_chain(input_file, ligand_chain_ids, verbose=True):
    """Maps original chain IDs to PDBFixer internal representation."""
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("orig", input_file)
    original_chain_ids = [chain.id for model in structure for chain in model]
    verbose and print(f"• Original chain IDs: {original_chain_ids}")

    fixer = PDBFixer(filename=input_file)
    chain_map = {orig_chain.id: new_chain for orig_chain, new_chain in zip(fixer.topology.chains(), original_chain_ids)}
    # verbose and print(f"• Map PDBFixer:OriginalPDB Chain Ids. {chain_map}")
    ligand_chain_ids = [chain_map[orig_id] for orig_id in ligand_chain_ids]
    
    verbose and print(f"• Ligand chain ids: {ligand_chain_ids}")
    return ligand_chain_ids

def get_ligand_atoms(fixer, ligand_chain_ids):
    """Gets atoms from specified ligand chains."""
    ligand_atoms = []
    for chain in fixer.topology.chains():
        if chain.id in ligand_chain_ids:
            for residue in chain.residues():
                for atom in residue.atoms():
                    ligand_atoms.append(atom)
    return ligand_atoms

def get_residues_near_ligand(fixer, ligand_atoms, distance_threshold):
    """Finds residues near ligand based on distance threshold."""
    near_residues = set()
    
    # Step 1: Convert ligand atom coordinates from nanometers to Ångströms
    ligand_coords = np.array([
        np.array(fixer.positions[atom.index].value_in_unit(unit.nanometer)) * 10.0
        for atom in ligand_atoms
    ])

    # Step 2: Iterate over all chains and residues to find those near ligand
    for chain in fixer.topology.chains():
        for residue in chain.residues():
            for atom in residue.atoms():
                atom_coord = np.array(fixer.positions[atom.index].value_in_unit(unit.nanometer)) * 10.0
                # Calculate distances between this atom and all ligand atoms
                distances = np.linalg.norm(ligand_coords - atom_coord, axis=1)
                if np.any(distances < distance_threshold):
                    near_residues.add(residue)
                    break  # Only need to find one atom within threshold
                    
    return near_residues

def get_chains_near_ligand(fixer, near_residues, ligand_chain_ids, verbose=True):
    """Gets chain IDs for residues near ligand, excluding ligand chains."""
    chain_ids = set()
    for chain in fixer.topology.chains():
        if chain.id in ligand_chain_ids:
            continue
        for residue in chain.residues():
            if residue in near_residues:
                chain_ids.add(chain.id)
                break  # Found at least one residue in this chain; move to the next chain.
    verbose and print(f"• Chains forming the target site: {sorted(chain_ids)}")  
    return sorted(set(chain_ids))