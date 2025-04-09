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
    from Bio.PDB.PDBParser import PDBParser
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
                record_type = line[0:6].strip()
                chain_search_key = (chain_id, record_type)
                if chain_search_key in target_chains:
                    fout.write(line)
            elif line.startswith("TER") and len(line) > 21:
                chain_id = line[21]
                if chain_id in target_chains:
                    fout.write(line)
            else:
                fout.write(line)
    
    rel_path = os.path.relpath(output_filepath)
    print(f"✅ Filtered PDB saved to {rel_path} containing chains {target_chains}")

def save_fixer_structure_to_pdb(pdb_structure, output_filepath, verbose=True):
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
    current_chain_ids = [chain.id for model in structure for chain in model]

    # If no target chains were specified, use all chains from the input file
    if target_chains is None:
        target_chains = current_chain_ids
    else:
        target_chains = [chain[0] for chain in target_chains]

    # Create a mapping between the chains in the input file and the target chains
    # This assumes chains are in the same order in both lists (PDBFixer typically preserves order)
    chain_mapping = {current_chain: correct_chain for current_chain, correct_chain in zip(current_chain_ids, target_chains)}
    verbose and print(f"• Chain mapping: {chain_mapping}")
    
    with open(input_filepath, 'r') as fin, open(output_filepath, 'w') as fout:
         for line in fin:
             # Process only ATOM and HETATM lines which contain atom coordinates
             if line.startswith("ATOM") or line.startswith("HETATM"):
                 # In PDB format, the chain ID is at position 21 (0-indexed)
                 chain_id = line[21:22]
                 if chain_id in chain_mapping:
                     # Get the target chain ID from the mapping
                     fixed_chain = chain_mapping[chain_id]
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

def save_clean_structure(input_file, output_file, skip_hetatm=True, verbose=True):
    """Loads PDB file, removes heteroatoms, and returns a PDBFixer object."""
    # Validate input file
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"File not found. {input_file}")
    if verbose:
        # Option 1: Get path relative to current working directory
        rel_path = os.path.relpath(input_file)
        print("\nPath to PDB file:", rel_path)
    
    try:
        # Step 1: Read PDB content
        with open(input_file, 'r') as f:
            lines = f.readlines()

        # Step 2: Process lines and filter HETATM records
        cleaned_lines = []
        chains = []
        chains_seen = set()
        # Skip water molecules only and keep other HETATM
        water_molecules = ["HOH", "WAT", "H2O", "SOL", "DOD"]

        for line in lines:
            if line.startswith("HETATM"):
                if skip_hetatm:  
                    # Skip all HETATM lines
                    continue
                else:
                    # If not skipping all HETATM, still skip water molecules
                    if any(water in line for water in water_molecules):
                        continue
            
            cleaned_lines.append(line)
            _count_chains(line, chains, chains_seen)
        
        # Step 3: Write cleaned content to file
        cleaned_pdb_text = "".join(cleaned_lines)
        with open(output_file, 'w') as f:
            f.write(cleaned_pdb_text)

        if verbose:
            rel_path = os.path.relpath(output_file)
            print(f"✅ Cleaned PDB saved to {rel_path} with chains: {chains}")
            
        return chains

    except Exception as e:
        print(f"❌ Error processing PDB file: {e}")
        raise

def _count_chains(line, chains, chains_seen):
    # Track chain IDs for atoms
    if line.startswith("ATOM") or line.startswith("HETATM"):
        chain_id = line[21:22]
        if chain_id not in chains_seen:
            chains.append(chain_id)
            chains_seen.add(chain_id)

    # Reset chains_seen when a TER record is encountered
    elif line.startswith("TER"):
        chains_seen.clear()

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
        if len(missing_residues_in_target_chains) > 0:
            print(f"\nMissing residues by target chains:")
            print(f"\n{'Chain ID':<10} {'Missing residues':<15} {'Position':<10}")
            print(f"{'-'*10} {'-'*15} {'-'*10}")
            for summary in summaries:
                print(f"{summary['chain_id']:<10} {summary['missing_residues_count']:<15} {summary['residue_index']:<10}")
        else:
            print("No missing residues found in target chains")

    return missing_residues_in_target_chains

#==============================================================================
#                      CHAIN AND RESIDUE ANALYSIS
#==============================================================================

def get_chain_indices(chain_map, chain_ids, record_type):
    if isinstance(chain_ids, str):
        result = [chain_ids]
    else:
        # Get all valid indices as a list
        result = [chain_map.get((chain_id, record_type)).get('fixer_index') 
                for chain_id in chain_ids 
                if (chain_id, record_type) in chain_map]
    print(f"• Chain indices: {result}")
    return result


def get_chains_to_extract(chain_map, chain_indices, fixer_object, cutoff, filter_type=None):
    selected_chains = None

    # Map the original chain IDs to the potentially new IDs after cleaning
    # chain_ids = map_orig_to_fixer_chain_ids(input_file, selection)
    # print(f"• Chain ids: {chain_ids}")
    
    # Identify atoms belonging to the specified HETATM chains
    ligand_atoms = get_ligand_atoms(fixer_object, chain_indices, filter_type) 
  
    # Find protein residues within the cutoff distance of the HETATM atoms
    residues_near_ligand = get_residues_near_ligand(
        fixer_object, ligand_atoms, distance_threshold=cutoff
    )

    print(f"• Residues near ligand: {len(residues_near_ligand)}")

    # Determine which chains contain the identified nearby residues
    selected_chain_indices = get_chains_near_ligand(fixer_object, residues_near_ligand, chain_indices)

    # Reverse map fixer chain index to fixer chain id. 
    selected_chain_ids = [(k[0], k[1]) for k, v in chain_map.items() if v.get('fixer_index') in selected_chain_indices]

    print(f"• List of chain indices near ligand: {sorted(selected_chain_ids)}")  

    return selected_chain_ids


def count_chains_biopython(input_file):
      # Use Biopython's PDB parser for validation
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("structure", input_file)
    # Validate against Biopython's model
    print("\nBiopython structure analysis:")
    for model in structure:
        for chain in model:
            chain_id = chain.id
            standard_residues = 0
            nonstandard_residues = 0
            
            for residue in chain:
                if residue.id[0] == ' ':
                    standard_residues += 1
                else:
                    nonstandard_residues += 1
                    
            print(f"Chain {chain_id}: {standard_residues} standard, {nonstandard_residues} non-standard residues")


def get_original_chain_ids(input_file):
    """Extract chain IDs directly from PDB file, preserving duplicates."""
    chain_ids = []
    seen_chains = set()
    
    # Track counts for debugging
    atom_chains = set()
    hetatm_chains = set()
    ter_count = 0
    
    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM"):
                chain_id = line[21:22].strip()  # Chain ID is at position 21 (0-indexed)
                atom_chains.add(chain_id)
                
                # Only add if it's a new chain (not seen since the last one)
                if chain_id not in seen_chains:
                    chain_ids.append(chain_id)
                    seen_chains.add(chain_id)
                
            elif line.startswith("HETATM"):
                chain_id = line[21:22].strip()
                hetatm_chains.add(chain_id)
                
                # Only add if it's a new chain (not seen since the last one)
                if chain_id not in seen_chains:
                    chain_ids.append(chain_id)
                    seen_chains.add(chain_id)
                    
            # Clear seen_chains when we encounter a TER record
            # This ensures that if the same chain ID appears again, after a TER record, it will be added to the list again
            elif line.startswith("TER"):
                ter_count += 1
                seen_chains.clear()
    
    # Print debugging information
    print("\n=== PDB Chain Analysis ===")
    print(f"ATOM chains: {sorted(atom_chains)}")
    print(f"HETATM chains: {sorted(hetatm_chains)}")
    print(f"TER records: {ter_count}")
    print(f"Chain IDs extracted (with duplicates from TER): {chain_ids}")
    print(f"Unique chain IDs: {sorted(set(chain_ids))}")
    
    return chain_ids


def get_fixer_chain_ids(input_file):
    """Extract chain IDs directly from PDBFixer object."""
    fixer = PDBFixer(filename=input_file)
    fixer_chain_ids = [chain.id for chain in fixer.topology.chains()]
    return fixer, fixer_chain_ids


def map_orig_to_fixer_chain_ids(input_file, ligand_chain_ids, verbose=True):
    """Maps original chain IDs to PDBFixer internal representation."""
    # parser = PDB.PDBParser(QUIET=True)
    # structure = parser.get_structure("orig", input_file)
    # original_chain_ids = [chain.id for model in structure for chain in model]
    original_chain_ids = get_original_chain_ids(input_file)
    verbose and print(f"• Original chain IDs: {original_chain_ids}")

    fixer, fixer_chain_ids = get_fixer_chain_ids(input_file)
    print(f"• Fixer chain ids: {fixer_chain_ids}")
    
    chain_map = {orig_chain.id: new_chain.index for orig_chain, new_chain in zip(fixer.topology.chains(), original_chain_ids)}
    print(f"• Chain map: {chain_map}")
    # verbose and print(f"• Map PDBFixer:OriginalPDB Chain Ids. {chain_map}")
    ligand_chain_ids = [chain_map[orig_id] for orig_id in ligand_chain_ids]
    
    verbose and print(f"• Ligand chain ids: {ligand_chain_ids}")
    return ligand_chain_ids

def get_ligand_atoms(fixer, chain_indices, filter_type=None):
    """Gets atoms from specified ligand chains."""

    ligand_atoms = []
    for chain in fixer.topology.chains():
        if chain.index in chain_indices:
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

def get_chains_near_ligand(fixer, near_residues, chain_indices, verbose=True):
    """Gets chain IDs for residues near ligand, excluding ligand chains."""
    chain_ids = set()
    for chain in fixer.topology.chains():
        if chain.index in chain_indices:
            continue
        for residue in chain.residues():
            if residue in near_residues:
                chain_ids.add(chain.index)
                break  # Found at least one residue in this chain; move to the next chain.
    
    return sorted(set(chain_ids))

def map_pdbfixer_chains_to_original(pdbfixer_obj, original_chain_ids):
    """
    Simply counts the number of chains in a PDBFixer object.
    """
    std_residues = set([
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 
        'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
    ])

    chain_count = 0
    chain_ids = {}
    for chain in pdbfixer_obj.topology.chains():
        
        if chain._residues[0].name in std_residues:
            chain_ids[(chain.id, 'ATOM')] = {'fixer_index':chain_count, 'original_chain_id': original_chain_ids[chain_count]}
        else:
            chain_ids[(chain.id, 'HETATM')] = {'fixer_index':chain_count, 'original_chain_id': original_chain_ids[chain_count]}
        chain_count += 1

    print(f"PDBFixer contains {chain_count} chains")
    return chain_ids