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
    python run.py --file_prefix <PDB_ID> --input_file <input.pdb> --include_chains A,B --verbose
    
    or
    
    python run.py --file_prefix <PDB_ID> --input_file <input.pdb> --reference_atom_chains C --verbose

Options:
    -file_prefix                    Prefix of the output files (required)
    --input_file                 Path to input PDB file (required)
    --include_chains             Comma-separated list of protein chains to be directly included in processing (e.g., A,B)
    --reference_atom_chains      Comma-separated list of ligand chains used as reference to identify nearby protein chains
    --reference_hetatm_chains    Comma-separated list of HETATM chains used as reference to identify nearby protein chains
    --cutoff                     Distance cutoff (Å) for identifying binding site residues (default: 5.0)
    --ph                         pH value for protonation (default: 7.4)
    --output_dir                 Directory for output files (default: results)
    --verbose                    Enable verbose output for detailed processing information

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
    parser.add_argument("--input_file", required=True, help="Path to input PDB file")
    
    # Chain selection (mutually exclusive group)
    chain_group = parser.add_argument_group("Chain Selection Options (specify one)")
    selection = chain_group.add_mutually_exclusive_group()
    selection.add_argument("--include_chains", help="Comma-separated list of protein chains to be directly included in processing (e.g., A,B)")
    selection.add_argument("--reference_atom_chains", help="Comma-separated list of ligand chains used as reference to identify nearby protein chains")
    selection.add_argument("--reference_hetatm_chains", help="Comma-separated list of HETATM chains used as reference to identify nearby protein chains")
    
    # Additional parameters
    parser.add_argument("--cutoff", type=float, default=5.0, help="Distance cutoff (Å) for identifying binding site residues")
    parser.add_argument("--ph", type=float, default=7.4, help="pH value for protonation")
    parser.add_argument("--output_dir", default="results", help="Directory for output files")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output for detailed processing information")
    parser.add_argument("--skip_molprobity", action="store_true", help="Disable molprobity and side chain optimization")
    parser.add_argument("--file_prefix", default="xx", help="Example PDB identifier (e.g., 1JY7)")
    
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
    VERBOSE = args.verbose
    FILE_PREFIX = args.file_prefix
    CUTOFF_ANGSTROM = args.cutoff
    RESULTS_FOLDER = args.output_dir
    SKIP_MOLPROBITY = args.skip_molprobity
    pH_VALUE = args.ph

    # Parse chain IDs from command line arguments
    INCLUDE_CHAIN_IDS = args.include_chains.split(',') if args.include_chains else None
    REFERENCE_ATOM_CHAIN_IDS = args.reference_atom_chains.split(',') if args.reference_atom_chains else None
    REFERENCE_HETATM_CHAIN_IDS = args.reference_hetatm_chains.split(',') if args.reference_hetatm_chains else None
    print(f"REFERENCE_ATOM_CHAIN_IDS: {REFERENCE_ATOM_CHAIN_IDS}")
    print(f"REFERENCE_HETATM_CHAIN_IDS: {REFERENCE_HETATM_CHAIN_IDS}")

    # Use a relative path based on the script location
    script_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    CONFIG_FILE = os.path.join(script_dir, 'scripts', 'config_env.json')
        
    # Store file paths for final summary
    file_paths = {}
    results_folder = os.path.abspath(RESULTS_FOLDER)
    if not os.path.exists(RESULTS_FOLDER):
        os.makedirs(RESULTS_FOLDER)
        VERBOSE and print(f"• Results folder created: {RESULTS_FOLDER}")
    else:
        VERBOSE and print(f"• Results folder exists: {RESULTS_FOLDER}")
        
    # Store all file paths
    INPUT_ORIGINAL_FILE = args.input_file
    file_paths['input'] = INPUT_ORIGINAL_FILE
    CLEANED_FILE = os.path.join(results_folder, f"{FILE_PREFIX}_0_structure_cleaned.pdb")
    file_paths['cleaned'] = CLEANED_FILE
    CLEANED_FILE2 = os.path.join(results_folder, f"{FILE_PREFIX}_0_structure_cleaned2.pdb")
    file_paths['cleane2d'] = CLEANED_FILE2
    SELECTED_CHAINS_FILE = os.path.join(results_folder, f"{FILE_PREFIX}_0_structure_selected_chains.pdb")
    file_paths['binding_site'] = SELECTED_CHAINS_FILE
    COMPLETED_TEMP_FILE = os.path.join(results_folder, f"{FILE_PREFIX}_1_structure_completed_temp.pdb")
    file_paths['temp_refined'] = COMPLETED_TEMP_FILE
    COMPLETED_FINAL_FILE = os.path.join(results_folder, f"{FILE_PREFIX}_2_structure_completed_final.pdb")
    file_paths['refined'] = COMPLETED_FINAL_FILE
    FLIPPED_H_TEMP_FILE = os.path.join(results_folder, f"{FILE_PREFIX}_3_structure_flipped_h_temp.pdb")
    file_paths['optimized'] = FLIPPED_H_TEMP_FILE
    FLIPPED_H_FINAL_FILE = os.path.join(results_folder, f"{FILE_PREFIX}_4_structure_flipped_h_final.pdb")
    file_paths['edited'] = FLIPPED_H_FINAL_FILE
    PROTONATED_PQR_FILE = os.path.join(results_folder, f"{FILE_PREFIX}_5_structure_protonated.pqr")
    file_paths['protonated_pqr'] = PROTONATED_PQR_FILE
    DOCKING_FILE = os.path.join(results_folder, f"{FILE_PREFIX}_6_structure_docking.pdbqt")
    file_paths['docking'] = DOCKING_FILE
    # OUTPUT_FINAL_PDB_FILE = os.path.join(results_folder, f"{PDB_ID}_6_structure_pka_protonated.pdb")
    # file_paths['final'] = OUTPUT_FINAL_PDB_FILE

    # Print initial message
    if not VERBOSE:
        print("\nProcessing structure preparation workflow...")
    
    #-----------------------------------------------------------------------
    # Step 0: PRINT SETTINGS
    #-----------------------------------------------------------------------
    add_separator("STEP 0: Settings")
    print(f"CONFIG_FILE: {CONFIG_FILE}")
    print(f"VERBOSE: {VERBOSE}")
    print(f"PDB_ID: {FILE_PREFIX}")
    print(f"INCLUDE_CHAIN_IDS: {INCLUDE_CHAIN_IDS}")
    print(f"REFERENCE_ATOM_CHAIN_IDS: {REFERENCE_ATOM_CHAIN_IDS}")
    print(f"REFERENCE_HETATM_CHAIN_IDS: {REFERENCE_HETATM_CHAIN_IDS}")
    print(f"CUTOFF_ANGSTROM: {CUTOFF_ANGSTROM}")
    print(f"pH_VALUE: {pH_VALUE}")
    print(f"RESULTS_FOLDER: {RESULTS_FOLDER}")
    print(f"INPUT_ORIGINAL_FILE: {INPUT_ORIGINAL_FILE}")
    

    #-----------------------------------------------------------------------
    # Step 1: Clean and return original structure as PDBFixer object
    #-----------------------------------------------------------------------
    VERBOSE and add_separator("STEP 1: LOADING AND CLEANING ORIGINAL STRUCTURE")
    skip_hetatm = True if REFERENCE_HETATM_CHAIN_IDS is None else False   
    print(f"• Skip HETATM: {skip_hetatm}")

    # Save clean structure and get original chain ids
    original_chain_ids = save_clean_structure(INPUT_ORIGINAL_FILE, CLEANED_FILE, skip_hetatm=skip_hetatm, verbose=VERBOSE)

    # Load structure as PDBFixer object
    fixer_original = load_structure_as_pdbfixer(CLEANED_FILE)

    # Map pdbfixer chain ids to chain indices and original chain ids
    chain_map = map_pdbfixer_chains_to_original(fixer_original, original_chain_ids)


    #-----------------------------------------------------------------------
    # Step 2: Define chains to be extracted
    #-----------------------------------------------------------------------
    VERBOSE and add_separator("STEP 2: DEFINING CHAINS TO BE EXTRACTED")
    # Three methods for chain selection, prioritized in this order:
    # 1. Direct inclusion: Explicitly include specified chains without proximity calculations
    # 2. Reference-based (HETATM): Use specified HETATM chains as reference points to find nearby protein chains
    # 3. Reference-based (atoms): Use specified ligand chains as reference points to find nearby protein chains

    chains_to_extract = None
    # Method 1: Directly include specified chains without proximity calculations
    if INCLUDE_CHAIN_IDS is not None:
        chains_to_extract = [(chain_id, "ATOM") for chain_id in INCLUDE_CHAIN_IDS]
    # Method 2: Use HETATM chains as reference points
    # (protein chains within cutoff distance of these reference HETATM chains will be selected)
    elif REFERENCE_HETATM_CHAIN_IDS is not None:
        chain_indices = get_chain_indices(chain_map, REFERENCE_HETATM_CHAIN_IDS, record_type="HETATM")
        chains_to_extract = get_chains_to_extract(chain_map, chain_indices, fixer_original, CUTOFF_ANGSTROM,  filter_type='HETATM')
        print(f"• Chains to extract NEAR HETATM CHAINS: {chains_to_extract}")
    # Method 3: Use Atom chains as reference points
    # (protein chains within cutoff distance of these reference ligand chains will be selected)
    elif REFERENCE_ATOM_CHAIN_IDS is not None:
        chain_indices = get_chain_indices(chain_map, REFERENCE_ATOM_CHAIN_IDS, record_type="ATOM")
        chains_to_extract = get_chains_to_extract(chain_map, chain_indices, fixer_original, CUTOFF_ANGSTROM,  filter_type='ATOM')
        print(f"• Chains to extract NEAR ATOM CHAINS: {chains_to_extract}")


    # #-----------------------------------------------------------------------
    # # Step 3: Extract chains
    # #-----------------------------------------------------------------------
    VERBOSE and add_separator("STEP 3: EXTRACTING CHAINS TO PDB")
    if chains_to_extract is not None:    
        extract_chains_to_pdb(CLEANED_FILE, SELECTED_CHAINS_FILE, chains_to_extract)
    else:
        # Default fallback: If no chains could be identified, use the entire structure
        SELECTED_CHAINS_FILE = CLEANED_FILE
        VERBOSE and print("No ligand chains specified, using the entire structure")
    

    #-----------------------------------------------------------------------
    # Step 4: Count missing residues
    #-----------------------------------------------------------------------
    VERBOSE and add_separator("STEP 4: COUNTING MISSING RESIDUES")
    original_pdbfixer = PDBFixer(filename=INPUT_ORIGINAL_FILE)
    original_pdbfixer.findMissingResidues()
    # We use the original structure to count missing residues, to make sure we are not missing any residues
    missing_residues_original_structure = get_missing_residues_by_chain(original_pdbfixer, chains_to_extract, verbose=VERBOSE)

    # From here on we use the selected chains. Clean structure again to exclude any remaining HETATM
    save_clean_structure(SELECTED_CHAINS_FILE, CLEANED_FILE, skip_hetatm=True, verbose=VERBOSE)
    

    if len(missing_residues_original_structure) > 0:
        #-----------------------------------------------------------------------
        # Step 5: Complete missing residues & atoms
        #-----------------------------------------------------------------------
        VERBOSE and add_separator("STEP 5: COMPLETE MISSING RESIDUES & ATOMS")

        pdb_fixer_object = load_structure_as_pdbfixer(CLEANED_FILE)
        
        completed_refined_fixer, residue_count_before, atom_count_before, residue_count_after, atom_count_after = complete_missing_structure(
            pdb_fixer_object, missing_residues_dict=missing_residues_original_structure, verbose=VERBOSE
        )

        #-----------------------------------------------------------------------
        # Step 6: SAVE FINAL (COMPLETE) STRUCTURE
        #-----------------------------------------------------------------------
        VERBOSE and add_separator("STEP 6: SAVE FINAL (COMPLETE) STRUCTURE")
        save_fixer_structure_to_pdb(completed_refined_fixer, COMPLETED_TEMP_FILE, verbose=VERBOSE)

        # PDBFixer renames chains. We restore it to the original chain IDs
        restore_original_chain_ids(COMPLETED_TEMP_FILE, COMPLETED_FINAL_FILE, chains_to_extract, verbose=VERBOSE)

        print(f"\nStructure completion statistics:")
        print(f"  {'Description':<30} {'Before':<10} {'After':<10}")
        print(f"  {'-'*30} {'-'*10} {'-'*10}")
        print(f"  {'Residue count':<30} {residue_count_before:<10} {residue_count_after:<10}")
        print(f"  {'Atom count':<30} {atom_count_before:<10} {atom_count_after:<10}")
        print(f"  {'Added atoms':<30} {'':<10} {atom_count_after - atom_count_before:<10}")

    else:
        VERBOSE and add_separator("STEP 5: SKIPPING STRUCTURE COMPLETION")
        restore_original_chain_ids(CLEANED_FILE, COMPLETED_FINAL_FILE, chains_to_extract, verbose=VERBOSE)

   

    if SKIP_MOLPROBITY:
        #-----------------------------------------------------------------------
        # Step 6: Optimize structure for docking
        #-----------------------------------------------------------------------
        # We skip molprobity and side chain optimization and instead use pdb2pqr to protonate the structure. 
        # Then we use MGLTools to create the PDBQT file.
        VERBOSE and add_separator("STEP 7: SKIPPING MOLPROBITY")
        VERBOSE and add_separator("STEP 8: PROTONATION WITH PDB2PQR")
        run_program("PDB2PQR", COMPLETED_FINAL_FILE, PROTONATED_PQR_FILE, verbose=VERBOSE, pH_value=pH_VALUE, config_file=CONFIG_FILE)

        #-----------------------------------------------------------------------
        # Step 7: Create PDBQT file for AutoDock Vina
        #-----------------------------------------------------------------------
        VERBOSE and add_separator("STEP 8: CREATING PDBQT FILE FOR AutoDock Vina")  
        run_program("MGLTools", PROTONATED_PQR_FILE, DOCKING_FILE, verbose=VERBOSE, pH_value=pH_VALUE, config_file=CONFIG_FILE)
    else:
        #-----------------------------------------------------------------------
        # Step 6: Optimize structure with MolProbity
        #-----------------------------------------------------------------------
        VERBOSE and add_separator("STEP 7: OPTIMIZING STRUCTURE FOR DOCKING")
        run_program("MolProbity", COMPLETED_FINAL_FILE, FLIPPED_H_TEMP_FILE, verbose=VERBOSE, pH_value=pH_VALUE, config_file=CONFIG_FILE)

        #-----------------------------------------------------------------------
        # Step 7: Convert Molprobity PDB file to proper PDB with OpenBabel
        #-----------------------------------------------------------------------
        VERBOSE and add_separator("STEP 5: POSTPROCESSING MOLPPROBITY FILE")   
        run_program("OpenBabel", FLIPPED_H_TEMP_FILE, FLIPPED_H_FINAL_FILE, verbose=VERBOSE, pH_value=pH_VALUE, config_file=CONFIG_FILE)

        #-----------------------------------------------------------------------
        # Step 8: Pronoate with pdb2pqr for correct protonation
        #-----------------------------------------------------------------------
        VERBOSE and add_separator("STEP 6: PROTONATION WITH PDB2PQR")
        run_program("PDB2PQR", FLIPPED_H_FINAL_FILE, PROTONATED_PQR_FILE, verbose=VERBOSE, pH_value=pH_VALUE, config_file=CONFIG_FILE)

        #-----------------------------------------------------------------------
        # Step 9: Create PDBQT file for AutoDock Vina
        #-----------------------------------------------------------------------
        VERBOSE and add_separator("STEP 7: CREATING PDBQT FILE FOR AutoDock Vina")  
        run_program("MGLTools", PROTONATED_PQR_FILE, DOCKING_FILE, verbose=VERBOSE, pH_value=pH_VALUE, config_file=CONFIG_FILE)

    
    # Print output files summary in non-verbose mode
    print("\nOutput files:")
    # Using string formatting with fixed width fields for consistent columns
    print(f"  {'File Type':<30} {'Path':<50}")
    print(f"  {'-'*30} {'-'*50}")
    print(f"  {'Binding site structure':<30} {os.path.relpath(SELECTED_CHAINS_FILE):<50}")
    print(f"  {'Refined structure':<30} {os.path.relpath(COMPLETED_FINAL_FILE):<50}")
    print(f"  {'Optimized structure':<30} {os.path.relpath(FLIPPED_H_FINAL_FILE):<50}")
    print(f"  {'Final PDBQT for docking':<30} {os.path.relpath(DOCKING_FILE):<50}")
    print("\n✅ Structure preparation completed successfully!")


if __name__ == "__main__":
    main()