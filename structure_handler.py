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
import argparse
from subprocess_handler import optimize_structure_with_molprobity, convert_to_pdbqt_with_mgltools, convert_file_with_openbabel

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

def load_clean_structure(pdb_filepath, no_hetatm=True, verbose=True):
    """Loads PDB file, removes heteroatoms, and returns a PDBFixer object."""
    # Validate input file
    if not os.path.exists(pdb_filepath):
        raise FileNotFoundError(f"File not found. {pdb_filepath}")
    if verbose:
        # Option 1: Get path relative to current working directory
        rel_path = os.path.relpath(pdb_filepath)
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
    
    # Prepare output filename
    base, ext = os.path.splitext(pdb_filepath)
    cleaned_filepath = f'{base}_clean{ext}'
    
    try:
        # Step 1: Read and filter PDB content
        with open(pdb_filepath, 'r') as f:
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
        with open(cleaned_filepath, 'w') as f:
            f.write(cleaned_pdb_text)
        
        # Step 4: Create and return PDBFixer object
        fixer = PDBFixer(filename=cleaned_filepath)
        if verbose:
            rel_cleaned_path = os.path.relpath(cleaned_filepath)
            print(f"✅ Cleaned structure from HETATM and loaded PDBFixer object from: {rel_cleaned_path}")
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

#==============================================================================
#                          MAIN EXECUTION
#==============================================================================

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Process and prepare protein structures for molecular docking.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required parameter
    parser.add_argument("--pdb_id", required=True, help="PDB identifier (e.g., 1JY7)")
    parser.add_argument("--input_file", required=True, help="Path to input PDB file")
    
    # Chain selection (mutually exclusive group)
    chain_group = parser.add_argument_group("Chain Selection Options (specify one)")
    selection = chain_group.add_mutually_exclusive_group()
    selection.add_argument("--target_chains", help="Comma-separated list of target protein chains (e.g., A,B)")
    selection.add_argument("--ligand_chains", help="Comma-separated list of ligand chains for identifying binding site")
    selection.add_argument("--hetatm_chains", help="Comma-separated list of HETATM chains for identifying binding site")
    
    # Additional parameters
    parser.add_argument("--cutoff", type=float, default=5.0, 
                        help="Distance cutoff (Å) for identifying binding site residues")
    parser.add_argument("--ph", type=float, default=7.4, 
                        help="pH value for protonation")
    parser.add_argument("--output_dir", default=".", 
                        help="Directory for output files")
    
    # Parse the arguments
    args = parser.parse_args()
    
    # Validate input file exists
    if not os.path.exists(args.input_file):
        parser.error(f"Input file not found: {args.input_file}")
    
    # # Validate that at least one chain selection method is provided
    # if not any([args.target_chains, args.ligand_chains, args.hetatm_chains]):
    #     parser.error("You must specify at least one of --target_chains, --ligand_chains, or --hetatm_chains")
    
    # Validate cutoff is positive
    if args.cutoff <= 0:
        parser.error("Distance cutoff must be a positive value")
    
    # Validate pH is in reasonable range
    if not (0 <= args.ph <= 14):
        parser.error("pH value must be between 0 and 14")
    
    return args

def main():
    """
    Execute the main workflow with validated arguments.
    
    This function implements the structure preparation workflow:
    1. Loads and cleans original structure
    2. Extracts binding site
    3. Refines structure by completing missing elements
    4. Optimizes structure for docking
    5. Creates PDBQT file for AutoDock Vina
    """
    try:
        # Add verbose flag to control output level
        VERBOSE = True # Set to False for minimal output, True for full detail
        PDB_ID = "2v7q"
        TARGET_CHAIN_IDS = "J"      # Seit either target chains or Ligand chains
        LIGAND_CHAIN_IDS = None  # Chain IDs that belong to the ligand
        HETATM_CHAIN_IDS = None # Chain IDs that belong to HETATM
        CUTOFF_ANGSTROM = 5.0     # Distance cutoff for identifying binding site residues
        pH_VALUE = 7.4            # pH for protonation
        
        # Store file paths for final summary
        file_paths = {}
        
        # Define file paths for the workflow
        RESULTS_FOLDER = f'results'

        results_folder = os.path.abspath(RESULTS_FOLDER)
        if not os.path.exists(RESULTS_FOLDER):
            os.makedirs(RESULTS_FOLDER)
            VERBOSE and print(f"• Results folder created: {RESULTS_FOLDER}")
        else:
            VERBOSE and print(f"• Results folder exists: {RESULTS_FOLDER}")
           
        # Store all file paths
        INPUT_ORIGINAL_FILE = os.path.join(results_folder, f'{PDB_ID}_original.pdb')
        file_paths['input'] = INPUT_ORIGINAL_FILE
        OUTPUT_BINDING_SITE_FILE = os.path.join(results_folder, f"{PDB_ID}_1_binding_site.pdb")
        file_paths['binding_site'] = OUTPUT_BINDING_SITE_FILE
        OUTPUT_TEMP_REFINED_FILE = os.path.join(results_folder, f"{PDB_ID}_2_binding_site_temp_refined.pdb")
        file_paths['temp_refined'] = OUTPUT_TEMP_REFINED_FILE
        OUTPUT_REFINED_FILE = os.path.join(results_folder, f"{PDB_ID}_3_binding_site_refined.pdb")
        file_paths['refined'] = OUTPUT_REFINED_FILE
        OUTPUT_FLIPPED_H_FILE = os.path.join(results_folder, f"{PDB_ID}_4_binding_site_refined_flipped_h.pdb")
        file_paths['optimized'] = OUTPUT_FLIPPED_H_FILE
        OUTPUT_FLIPPED_H_EDITED_FILE = os.path.join(results_folder, f"{PDB_ID}_5_binding_site_refined_flipped_h_edited.pdb")
        file_paths['edited'] = OUTPUT_FLIPPED_H_EDITED_FILE
        OUTPUT_PDBQT_DOCKING_FILE = os.path.join(results_folder, f"{PDB_ID}_6_binding_site_refined_docking.pdbqt")
        file_paths['docking'] = OUTPUT_PDBQT_DOCKING_FILE
        OUTPUT_PROTONATED_PQR_FILE = os.path.join(results_folder, f"{PDB_ID}_7_binding_site_protonated.pqr")
        file_paths['protonated_pqr'] = OUTPUT_PROTONATED_PQR_FILE
        OUTPUT_FINAL_PDB_FILE = os.path.join(results_folder, f"{PDB_ID}_8_binding_site_protonated.pdb")
        file_paths['final'] = OUTPUT_FINAL_PDB_FILE
 
        # Print initial message
        if not VERBOSE:
            print("\nProcessing structure preparation workflow...")
        
        #-----------------------------------------------------------------------
        # Step 1: Load and clean original structure
        #-----------------------------------------------------------------------
        VERBOSE and add_separator("STEP 1: LOADING AND CLEANING ORIGINAL STRUCTURE")
        no_hetatm = True if HETATM_CHAIN_IDS is not None else False
        fixer_original = load_clean_structure(INPUT_ORIGINAL_FILE, no_hetatm=no_hetatm, verbose=VERBOSE)
        
        #-----------------------------------------------------------------------
        # Step 2: Extract binding site
        #-----------------------------------------------------------------------
        VERBOSE and add_separator("STEP 2: EXTRACTING BINDING SITE")
        target_chains = None
        if TARGET_CHAIN_IDS is not None:
            target_chains = TARGET_CHAIN_IDS
        elif HETATM_CHAIN_IDS is not None:
            hetatm_chain_ids = map_orig_to_new_chain(INPUT_ORIGINAL_FILE, HETATM_CHAIN_IDS)
            ligand_atoms = get_ligand_atoms(fixer_original, hetatm_chain_ids) 
            residues_near_ligand = get_residues_near_ligand(
                fixer_original, ligand_atoms, distance_threshold=CUTOFF_ANGSTROM
            )
            target_chains = get_chains_near_ligand(fixer_original, residues_near_ligand, hetatm_chain_ids)
        elif LIGAND_CHAIN_IDS is not None:
            ligand_chain_ids = map_orig_to_new_chain(INPUT_ORIGINAL_FILE, LIGAND_CHAIN_IDS)
            ligand_atoms = get_ligand_atoms(fixer_original, ligand_chain_ids) 
            residues_near_ligand = get_residues_near_ligand(
                fixer_original, ligand_atoms, distance_threshold=CUTOFF_ANGSTROM
            )
            target_chains = get_chains_near_ligand(fixer_original, residues_near_ligand, ligand_chain_ids)
                
        if target_chains is not None:    
            extract_chains_to_pdb(INPUT_ORIGINAL_FILE, OUTPUT_BINDING_SITE_FILE, target_chains)
        else:
            OUTPUT_BINDING_SITE_FILE = INPUT_ORIGINAL_FILE
            VERBOSE and print("No ligand chains specified, using the entire structure")
        
        #-----------------------------------------------------------------------
        # Step 3: Count missing residues
        #-----------------------------------------------------------------------
        VERBOSE and add_separator("STEP 3: REFINING STRUCTURE")
        original_pdbfixer = PDBFixer(filename=INPUT_ORIGINAL_FILE)
        original_pdbfixer.findMissingResidues()
        missing_residues_original_structure = get_missing_residues_by_chain(original_pdbfixer, target_chains, verbose=VERBOSE)

        #-----------------------------------------------------------------------
        # Step 4: Refine structure
        #-----------------------------------------------------------------------
        fixer_subset = load_clean_structure(OUTPUT_BINDING_SITE_FILE, verbose=VERBOSE)
        completed_refined_fixer, residue_count_before, atom_count_before, residue_count_after, atom_count_after = complete_missing_structure(
            fixer_subset, missing_residues_dict=missing_residues_original_structure, verbose=VERBOSE
        )
        save_structure_to_pdb(completed_refined_fixer, OUTPUT_TEMP_REFINED_FILE, verbose=VERBOSE)
        restore_original_chain_ids(OUTPUT_TEMP_REFINED_FILE, OUTPUT_REFINED_FILE, target_chains, verbose=VERBOSE)
        
        # Always print the residue counts, even in non-verbose mode
        if not VERBOSE:
            print(f"\nStructure completion statistics:")
            # Using string formatting with fixed width fields
            print(f"  {'Description':<30} {'Before':<10} {'After':<10}")
            print(f"  {'-'*30} {'-'*10} {'-'*10}")
            print(f"  {'Residue count':<30} {residue_count_before:<10} {residue_count_after:<10}")
            print(f"  {'Atom count':<30} {atom_count_before:<10} {atom_count_after:<10}")
            print(f"  {'Added atoms':<30} {'':<10} {atom_count_after - atom_count_before:<10}")

        #-----------------------------------------------------------------------
        # Step 4: Optimize structure for docking
        #-----------------------------------------------------------------------
        VERBOSE and add_separator("STEP 4: OPTIMIZING STRUCTURE FOR DOCKING")
        optimize_structure_with_molprobity(OUTPUT_REFINED_FILE, OUTPUT_FLIPPED_H_FILE, verbose=VERBOSE)

        #-----------------------------------------------------------------------
        # Step 5: Clean Molprobity file
        #-----------------------------------------------------------------------
        VERBOSE and add_separator("STEP 5: POSTPROCESSING MOLPPROBITY FILE")   
        convert_file_with_openbabel(OUTPUT_FLIPPED_H_FILE, OUTPUT_FLIPPED_H_EDITED_FILE, verbose=VERBOSE)

        #-----------------------------------------------------------------------
        # Step 6: Create PDBQT file for AutoDock Vina
        #-----------------------------------------------------------------------
        VERBOSE and add_separator("STEP 6: CREATING PDBQT FILE FOR AutoDock Vina")  
        convert_to_pdbqt_with_mgltools(OUTPUT_FLIPPED_H_EDITED_FILE, OUTPUT_PDBQT_DOCKING_FILE, verbose=VERBOSE)
        
        # Print output files summary in non-verbose mode
        if not VERBOSE:
            print("\nOutput files:")
            # Using string formatting with fixed width fields for consistent columns
            print(f"  {'File Type':<30} {'Path':<50}")
            print(f"  {'-'*30} {'-'*50}")
            print(f"  {'Binding site structure':<30} {os.path.relpath(OUTPUT_BINDING_SITE_FILE):<50}")
            print(f"  {'Refined structure':<30} {os.path.relpath(OUTPUT_REFINED_FILE):<50}")
            print(f"  {'Optimized structure':<30} {os.path.relpath(OUTPUT_FLIPPED_H_FILE):<50}")
            print(f"  {'Final PDBQT for docking':<30} {os.path.relpath(OUTPUT_PDBQT_DOCKING_FILE):<50}")
            print("\n✅ Structure preparation completed successfully!")
        
        # Commented out protonation step
        # #-------------------------------------------------------------------
        # # Step 6: Protonation for MD simulations
        # #-------------------------------------------------------------------
        # add_separator("STEP 6: PROTONATION FOR MD SIMULATIONS")
        # protonate_structure_with_pdb2pqr(output_flipped_h_edited_file, output_protonated_pqr_file, pH=ph_value)
        # convert_pqr_to_pdb(output_protonated_pqr_file, output_final_pdb_file)
    except Exception as e:
        print(f"\n❌ Error during structure preparation: {e}")

if __name__ == "__main__":
    main()