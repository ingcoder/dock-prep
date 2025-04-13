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

import sys
import traceback

try:
    import numpy as np
    from pdbfixer import PDBFixer

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
#                      CHAIN ID MAPPING
#==============================================================================

def map_pdbfixer_chains_to_original(pdbfixer_obj, original_chain_ids):
    """Maps PDBFixer chain IDs to original chain IDs."""
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
    print(f"✅ Created PDBFixer object with {chain_count} chains")
    
    return chain_ids


#==============================================================================
#                        CHAIN SELECTION
#==============================================================================

def get_chains_to_extract(chain_map, chain_ids, record_type, fixer_object, cutoff):
    selected_chains = None
    # Map the original chain IDs to the potentially new IDs after cleaning
    # chain_ids = map_orig_to_fixer_chain_ids(input_file, selection)
    # print(f"• Chain ids: {chain_ids}")

    chain_indices = _get_chain_indices(chain_map, chain_ids, record_type)
    
    # Identify atoms belonging to the specified reference molecule
    ligand_atoms = _get_ligand_atoms(fixer_object, chain_indices) 
    # Find protein residues within the cutoff distance of the reference molecule
    residues_near_ligand = _get_residues_near_ligand(fixer_object, ligand_atoms, distance_threshold=cutoff)
    print(f"• Residues near ligand: {len(residues_near_ligand)}")
    # Determine which chains contain the identified nearby residues
    selected_chain_indices = _get_chains_near_ligand(fixer_object, residues_near_ligand, chain_indices)
    # Reverse map fixer chain index to fixer chain id. 
    selected_chain_ids = [(k[0], k[1]) for k, v in chain_map.items() if v.get('fixer_index') in selected_chain_indices]
    return selected_chain_ids

def _get_chain_indices(chain_map, chain_ids, record_type):
    """Gets chain indices for specified chain IDs."""
    if isinstance(chain_ids, str):
        result = [chain_ids]
    else:
        # Get all valid indices as a list
        result = [chain_map.get((chain_id, record_type)).get('fixer_index') 
                for chain_id in chain_ids 
                if (chain_id, record_type) in chain_map]
    print(f"• Chain indices: {result}")
    return result

def _get_ligand_atoms(fixer, ligand_selection):
    """Gets atoms from specified ligand chains."""
    ligand_atoms = []
    
    # Case 1: ligand_selection is a specific residue (chain_index, residue_id)
    if isinstance(ligand_selection, tuple) and len(ligand_selection) == 2:
        chain_index, residue_id = ligand_selection
        for chain in fixer.topology.chains():
            if chain.index == chain_index:
                for residue in chain.residues():
                    if residue.id == residue_id:
                        for atom in residue.atoms():
                            ligand_atoms.append(atom)
    
    # Case 2: ligand_selection is a list of chain indices
    else:
        for chain in fixer.topology.chains():
            if chain.index in ligand_selection:
                for residue in chain.residues():
                    for atom in residue.atoms():
                        ligand_atoms.append(atom)
    
    return ligand_atoms

def _get_residues_near_ligand(fixer, ligand_atoms, distance_threshold):
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

def _get_chains_near_ligand(fixer, near_residues, chain_indices, verbose=True):
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


#==============================================================================
#                     HELPER FUNCTIONS
#==============================================================================

def load_structure_as_pdbfixer(input_file):
    """Creates and returns a PDBFixer object from a PDB file."""
    try:
        fixer = PDBFixer(filename=input_file)
        return fixer
    except Exception as e:
        print(f"Error: Failed to process PDB file with PDBFixer: {e}")
        sys.exit(1)

def count_residues(fixer, quiet=True, verbose=True):
    """Counts residues and atoms in a PDBFixer object."""
    residue_count = 0
    atom_count = 0
    for chain in fixer.topology.chains():
        for residue in chain.residues():
            residue_count += 1
            for atom in residue.atoms():
                atom_count += 1
    if quiet:
        return residue_count, atom_count
    else:
        print("• Residue count:", residue_count)
        print("• Atom count:", atom_count)
        return residue_count, atom_count
    
def count_chains_biopython(input_file):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("structure", input_file)

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


#==============================================================================
#                     MISSING RESIDUES ANALYSIS
#==============================================================================

def get_missing_residues_by_chain(original_pdbfixer, selected_chains, include_n_terminal_gaps=True, verbose=True):
    """Gets missing residues in specified chains."""
    summaries = []
    missing_residues_in_target_chains = {}

    selected_chains_ids = [chain[0] for chain in selected_chains] if selected_chains is not None else None

    original_pdbfixer.findMissingResidues()
    
    chain_index_to_object = {}
    for chain in original_pdbfixer.topology.chains():
        if selected_chains_ids is not None:
            if chain.id in selected_chains_ids:
                print(f"• Chain {chain.id} is in selected chains {selected_chains_ids}")
                chain_index_to_object[chain.index] = chain
        else:
            chain_index_to_object[chain.index] = chain
 
    chain_id_to_index = {}                
    # Step 2: Collect missing residues for target chains
    for (chain_index, residue_index), residues in original_pdbfixer.missingResidues.items():
        # Skip if chain not in our target chains
        if chain_index not in chain_index_to_object:
            continue
            
        # Skip residue index 0 (N-terminal gaps) unless explicitly requested
        if residue_index == 0 and not include_n_terminal_gaps:
            continue
            
        # Get chain and create summary
        chain = chain_index_to_object[chain_index]
        # verbose and print(f"Chain {chain_index}, residue_index: {residue_index}, residues: {residues}")

        chain_id_to_index[chain.id] = chain.index
        
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
        if summaries:
            print(f"\nMissing residues in selected chains:")
            print(f"\n{'Chain ID':<10} {'Missing residues':<15} {'Position':<10}")
            print(f"{'-'*10} {'-'*15} {'-'*10}")
            for summary in summaries:
                print(f"{summary['chain_id']:<10} {summary['missing_residues_count']:<15} {summary['residue_index']:<10}")
        else:
            print("No missing residues found in target chains")

    return missing_residues_in_target_chains, chain_id_to_index


#==============================================================================
#                        STRUCTURE REFINEMENT
#==============================================================================

def complete_missing_structure(pdb_fixer, found_missing_residues=None, found_missing_chain_id_to_index=None, verbose=True):
    """Completes missing residues and atoms in a PDBFixer object."""
    residue_count_before, atom_count_before = count_residues(pdb_fixer, quiet=True, verbose=verbose)
   
    missing_residues_dict = _transform_missing_residues(pdb_fixer, found_missing_residues, found_missing_chain_id_to_index)


    pdb_fixer.findMissingResidues()
    
    # Step 2: Handle missing residues (automatic or from dictionary)
    if missing_residues_dict is None:
        auto_missing_residues = pdb_fixer.missingResidues
        print(f"\n> Automatically detected missing residues: {len(auto_missing_residues)}")
        for chain in auto_missing_residues:
            print(f"• Chain {chain}: {len(auto_missing_residues[chain])} missing residue segments")
    else:        
        pdb_fixer.missingResidues = missing_residues_dict
    
    # Step 3: Identify and report missing atoms and terminals
    try:        # Find missing atoms (and internal structures, including terminals)
        pdb_fixer.findMissingAtoms()
        missing_atoms = pdb_fixer.missingAtoms
        missing_terminals = pdb_fixer.missingTerminals
        if verbose:
            print(f"> Missing atoms found in {len(missing_atoms)} residues")
            print(f"> Missing terminals found in {len(missing_terminals)} chains")
    except Exception as e:
        print(f"Error: Failed to find missing atoms: {e}")
        traceback.print_exc()
        sys.exit(1)
    
    # Step 4: Add missing atoms AND missing residues to the structure
    try:
        verbose and print("\n> Adding missing residues and atoms...")
        
        # This single call adds both missing residues AND missing atoms
        pdb_fixer.addMissingAtoms()
        
        # DEBUGGING: Check if residues were actually added
        topology_after = pdb_fixer.topology
        chain_count_after = sum(1 for _ in topology_after.chains())
        residue_count_check_after = sum(1 for _ in topology_after.residues())
    except Exception as e:
        print(f"Error: Failed to add missing atoms and residues: {e}")
        traceback.print_exc()
        sys.exit(1)
    # Hydrogens can be added here if needed:
    # pdb_fixer.addMissingHydrogens(pH=6.2)
    if verbose:
        print("\n> After structure completion:")
    
    residue_count_after, atom_count_after = count_residues(pdb_fixer, verbose=verbose)

    # Report changes
    if verbose:
        res_diff = residue_count_after - residue_count_before
        atom_diff = atom_count_after - atom_count_before
        print(f"> Added {res_diff} residues and {atom_diff} atoms to the structure")
        
    return pdb_fixer


def _transform_missing_residues(clean_fixer_object, found_missing_residues, found_missing_chain_id_to_index):
    """
    Translates missing residue information between different chain numbering systems.
    The function converts from chain indices to chain IDs and then back to chain indices in the new structure.
    It's properly taking the missing residues information from the original structure and transforming it to be compatible with the cleaned structure.
    
    Problem:
    When we fix the protein structure with addMissingResidues, the chains get renumbered.
    This means the chain indices we found earlier no longer match the new structure.
    
    Solution:
    This function creates a mapping between the original chain IDs and the new chain indices,
    allowing us to correctly apply the missing residue information to the right chains.
    
    Parameters:
    - pdb_fixer_object: The PDBFixer object with the cleaned structure
    - found_missing_residues: Dictionary of missing residues using original chain indices
    - found_missing_chain_id_to_index: Returned together with found_missing_residues and maps chain ids to indices
    
    Returns:
    Dictionary of missing residues with updated chain indices matching the cleaned structure
    """
    # Step 1. Maps chain ids to indices in the cleaned structure {B : 0, C : 1, ...}
    # ID -> Index
    cleaned_chain_id_to_index = {}
    for chain in clean_fixer_object.topology.chains():
        cleaned_chain_id_to_index[chain.id] = chain.index

    # Corrected missing residues dictionary
    transformed_missing_residues = {}

    # Get chain indices found in missing residues of original structure
    for (chain_idx, res_idx), residues in found_missing_residues.items():

        # Map indices to ids
        chain_id = None
        for orig_id, orig_idx in found_missing_chain_id_to_index.items():
            if orig_idx == chain_idx:
                chain_id = orig_id
                break
        
        # If id in cleaned structure map id to index of 
        if chain_id is not None and chain_id in cleaned_chain_id_to_index:
            new_chain_idx = cleaned_chain_id_to_index[chain_id]
            transformed_missing_residues[(new_chain_idx, res_idx)] = residues
            # if params['verbose']:
            #     print(f"• Mapped missing residues from chain index {chain_idx} to {new_chain_idx} (Chain ID: {chain_id})")

    return transformed_missing_residues
        
