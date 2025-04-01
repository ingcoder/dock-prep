"""
PDB to PDBQT Converter: Molecular Structure Preparation Pipeline

This script automates the preparation of protein structures for molecular docking and simulation.
It implements a complete workflow that processes PDB structures through several refinement stages:

1. Structure Loading and Cleaning: Reads a PDB file and performs initial cleanup
2. Binding Site Extraction: Identifies the binding site based on ligand proximity or specified chains
3. Structure Refinement: Completes missing residues and adds missing atoms
4. Structure Optimization: Optimizes hydrogen bonds and corrects histidine tautomers
5. PDBQT Conversion: Generates AutoDock Vina compatible PDBQT files with proper atom types and charges

The script can be used as a standalone tool with command-line arguments or programmatically
by modifying the parameters in the main() function.

Dependencies:
- PDBFixer for structure completion
- MolProbity for hydrogen optimization
- OpenBabel for structure cleaning
- MGLTools for PDBQT conversion

Usage:
    python run.py --pdb_id <PDB_ID> --input_file <input.pdb> --target_chains A,B
    
    or
    
    python run.py --pdb_id <PDB_ID> --input_file <input.pdb> --ligand_chains X --verbose

Options:
    --pdb_id            PDB identifier (required)
    --input_file        Path to input PDB file (required)
    --target_chains     Comma-separated list of target protein chains (e.g., A,B)
    --ligand_chains     Comma-separated list of ligand chains for identifying binding site
    --hetatm_chains     Comma-separated list of HETATM chains for identifying binding site
    --cutoff            Distance cutoff (Å) for identifying binding site residues (default: 5.0)
    --ph                pH value for protonation (default: 7.4)
    --output_dir        Directory for output files (default: results)
    --verbose           Enable verbose output for detailed processing information

Output:
    Multiple files representing different stages of the preparation process, including
    the final PDBQT file ready for molecular docking with AutoDock Vina.
"""

#==============================================================================
#                          MAIN EXECUTION
#==============================================================================

import argparse
import os
from dock_prep.subprocess_handler import run_program
from dock_prep.structure_handler import *

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
    parser.add_argument("--cutoff", type=float, default=5.0, help="Distance cutoff (Å) for identifying binding site residues")
    parser.add_argument("--ph", type=float, default=7.4, help="pH value for protonation")
    parser.add_argument("--output_dir", default="results", help="Directory for output files")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output for detailed processing information")
    
    args = parser.parse_args()
    
    # Validate input file exists
    if not os.path.exists(  os.path.abspath(args.input_file)):
        parser.error(f"Input file not found: {  os.path.abspath(args.input_file)}")
     
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
    
    # Parse command line arguments
    args = parse_arguments()
    CONFIG_FILE = os.path.join('/Users/ingrid/Projects/PdbqtConverter/dock-prep/scripts', 'config_env.json')    
    
    
    # Add verbose flag to control output level
    VERBOSE = args.verbose
    print(f"VERBOSE: {VERBOSE}")
    PDB_ID = args.pdb_id
    
    # Parse chain IDs from command line arguments
    TARGET_CHAIN_IDS = args.target_chains.split(',') if args.target_chains else None
    LIGAND_CHAIN_IDS = args.ligand_chains.split(',') if args.ligand_chains else None
    HETATM_CHAIN_IDS = args.hetatm_chains.split(',') if args.hetatm_chains else None
    
    # Get other parameters from arguments
    CUTOFF_ANGSTROM = args.cutoff
    pH_VALUE = args.ph
    
    # Store file paths for final summary
    file_paths = {}
    
    # Define file paths for the workflow
    RESULTS_FOLDER = args.output_dir

    results_folder = os.path.abspath(RESULTS_FOLDER)
    if not os.path.exists(RESULTS_FOLDER):
        os.makedirs(RESULTS_FOLDER)
        VERBOSE and print(f"• Results folder created: {RESULTS_FOLDER}")
    else:
        VERBOSE and print(f"• Results folder exists: {RESULTS_FOLDER}")
        
    # Store all file paths
    INPUT_ORIGINAL_FILE = args.input_file
    file_paths['input'] = INPUT_ORIGINAL_FILE
    OUTPUT_BINDING_SITE_FILE = os.path.join(results_folder, f"{PDB_ID}_0_binding_site.pdb")
    file_paths['binding_site'] = OUTPUT_BINDING_SITE_FILE
    OUTPUT_TEMP_REFINED_FILE = os.path.join(results_folder, f"{PDB_ID}_1_binding_site_temp_refined.pdb")
    file_paths['temp_refined'] = OUTPUT_TEMP_REFINED_FILE
    OUTPUT_REFINED_FILE = os.path.join(results_folder, f"{PDB_ID}_2_binding_site_refined.pdb")
    file_paths['refined'] = OUTPUT_REFINED_FILE
    OUTPUT_FLIPPED_H_FILE = os.path.join(results_folder, f"{PDB_ID}_3_binding_site_refined_flipped_h.pdb")
    file_paths['optimized'] = OUTPUT_FLIPPED_H_FILE
    OUTPUT_FLIPPED_H_EDITED_FILE = os.path.join(results_folder, f"{PDB_ID}_4_binding_site_refined_flipped_h_edited.pdb")
    file_paths['edited'] = OUTPUT_FLIPPED_H_EDITED_FILE
    OUTPUT_PDBQT_DOCKING_FILE = os.path.join(results_folder, f"{PDB_ID}_5_binding_site_refined_docking.pdbqt")
    file_paths['docking'] = OUTPUT_PDBQT_DOCKING_FILE
    OUTPUT_PROTONATED_PQR_FILE = os.path.join(results_folder, f"{PDB_ID}_6_binding_site_refined_protonated.pqr")
    file_paths['protonated_pqr'] = OUTPUT_PROTONATED_PQR_FILE
    OUTPUT_FINAL_PDB_FILE = os.path.join(results_folder, f"{PDB_ID}_7_binding_site_protonated.pdb")
    file_paths['final'] = OUTPUT_FINAL_PDB_FILE

    # Print initial message
    if not VERBOSE:
        print("\nProcessing structure preparation workflow...")
    
    #-----------------------------------------------------------------------
    # Step 0: Load and clean original structure
    #-----------------------------------------------------------------------
    add_separator("STEP 0: Settings")
    print(f"CONFIG_FILE: {CONFIG_FILE}")
    print(f"VERBOSE: {VERBOSE}")
    print(f"PDB_ID: {PDB_ID}")
    print(f"TARGET_CHAIN_IDS: {TARGET_CHAIN_IDS}")
    print(f"LIGAND_CHAIN_IDS: {LIGAND_CHAIN_IDS}")
    print(f"HETATM_CHAIN_IDS: {HETATM_CHAIN_IDS}")
    print(f"CUTOFF_ANGSTROM: {CUTOFF_ANGSTROM}")
    print(f"pH_VALUE: {pH_VALUE}")
    print(f"RESULTS_FOLDER: {RESULTS_FOLDER}")
    print(f"INPUT_ORIGINAL_FILE: {INPUT_ORIGINAL_FILE}")
    

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

    
    SKIP_MOLPROBITY = False
    if SKIP_MOLPROBITY:
        #-----------------------------------------------------------------------
        # Step 4: Optimize structure for docking
        #-----------------------------------------------------------------------
        VERBOSE and add_separator("STEP 6: SKIPPING MOLPROBITY")
        VERBOSE and add_separator("STEP 6: PROTONATION WITH PDB2PQR")
        run_program("PDB2PQR", OUTPUT_REFINED_FILE, OUTPUT_PROTONATED_PQR_FILE, verbose=VERBOSE, pH_value=pH_VALUE, config_file=CONFIG_FILE)

        #-----------------------------------------------------------------------
        # Step 6: Create PDBQT file for AutoDock Vina
        #-----------------------------------------------------------------------
        VERBOSE and add_separator("STEP 6: CREATING PDBQT FILE FOR AutoDock Vina")  
        run_program("MGLTools", OUTPUT_PROTONATED_PQR_FILE, OUTPUT_PDBQT_DOCKING_FILE, verbose=VERBOSE, pH_value=pH_VALUE, config_file=CONFIG_FILE)
    else:
        #-----------------------------------------------------------------------
        # Step 4: Optimize structure with MolProbity
        #-----------------------------------------------------------------------
        VERBOSE and add_separator("STEP 4: OPTIMIZING STRUCTURE FOR DOCKING")
        run_program("MolProbity", OUTPUT_REFINED_FILE, OUTPUT_FLIPPED_H_FILE, verbose=VERBOSE, pH_value=pH_VALUE, config_file=CONFIG_FILE)

        #-----------------------------------------------------------------------
        # Step 5: Convert Molprobity file to PDB with OpenBabel
        #-----------------------------------------------------------------------
        VERBOSE and add_separator("STEP 5: POSTPROCESSING MOLPPROBITY FILE")   
        run_program("OpenBabel", OUTPUT_FLIPPED_H_FILE, OUTPUT_FLIPPED_H_EDITED_FILE, verbose=VERBOSE, pH_value=pH_VALUE, config_file=CONFIG_FILE)

        #-----------------------------------------------------------------------
        # Step 6: Pronoate with pdb2pqr for correct protonation
        #-----------------------------------------------------------------------
        VERBOSE and add_separator("STEP 6: PROTONATION WITH PDB2PQR")
        run_program("PDB2PQR", OUTPUT_FLIPPED_H_EDITED_FILE, OUTPUT_PROTONATED_PQR_FILE, verbose=VERBOSE, pH_value=pH_VALUE, config_file=CONFIG_FILE)

        #-----------------------------------------------------------------------
        # Step 7: Create PDBQT file for AutoDock Vina
        #-----------------------------------------------------------------------
        VERBOSE and add_separator("STEP 6: CREATING PDBQT FILE FOR AutoDock Vina")  
        run_program("MGLTools", OUTPUT_PROTONATED_PQR_FILE, OUTPUT_PDBQT_DOCKING_FILE, verbose=VERBOSE, pH_value=pH_VALUE, config_file=CONFIG_FILE)

    
    
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
    
   

if __name__ == "__main__":
    main()