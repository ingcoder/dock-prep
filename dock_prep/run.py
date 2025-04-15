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
import os
from dock_prep.subprocess_handler import run_program
from dock_prep.structure_handler import *
from dock_prep.argument_handler import parse_arguments, extract_params_from_args
from dock_prep.file_handler2 import *

def run_pipeline(args):
    """Main workflow for preparing protein structures for docking."""
    #=======================================================================
    # STEP 0: INITIALIZATION
    #=======================================================================
    
    add_separator("STEP 0: INITIALIZATION")
    # Extract parameters from arguments
    params = extract_params_from_args(args)
    
    # Setup file paths
    file_paths, results_folder = setup_file_paths(
        params['file_prefix'], 
        params['results_folder'], 
        params['input_file']
    )
    
    # Print initial message
    if not params['verbose']:
        print("\nProcessing structure preparation workflow...")
    

    #=======================================================================
    # STEP 1: STRUCTURE CLEANING AND PREPARATION
    #=======================================================================
    
    params['verbose'] and add_separator("STEP 1: LOADING AND CLEANING ORIGINAL STRUCTURE")
    
    # Clean the structure and analyze chains
    original_chain_ids = save_clean_structure(
        file_paths['input'], 
        file_paths['cleaned'], 
        skip_hetatm=params['skip_hetatm'], 
        verbose=params['verbose']
    )
    
    # Load structure as PDBFixer object
    fixer_original = load_structure_as_pdbfixer(file_paths['cleaned'])

    # Map pdbfixer chain ids to original chain ids
    chain_map = map_pdbfixer_chains_to_original(fixer_original, original_chain_ids)


    #=======================================================================
    # STEP 2: CHAIN SELECTION
    #=======================================================================
    
    params['verbose'] and add_separator("STEP 2: DEFINING CHAINS TO BE EXTRACTED")
    # Three methods for chain selection, prioritized in this order:
    # 1. Direct inclusion: Explicitly include specified chains without proximity calculations
    # 2. Reference-based (HETATM): Use HETATM chains as reference points to find nearby protein chains
    # 3. Reference-based (atoms): Use peptide chains as reference points to find nearby protein chains

    chains_to_extract = None
    
    # Method 1: Direct inclusion of specified chains
    if params['include_chain_ids'] is not None:
        chains_to_extract = [(chain_id, "ATOM") for chain_id in params['include_chain_ids']]
    
    # Method 2 & 3: Reference-based selection using chain IDs
    elif params['reference_hetatm_chain_ids'] is not None or params['reference_atom_chain_ids'] is not None:
        # Determine which reference to use and its type
        if params['reference_hetatm_chain_ids'] is not None:
            reference_chain_ids = params['reference_hetatm_chain_ids']
            record_type = "HETATM"
        else:
            reference_chain_ids = params['reference_atom_chain_ids']
            record_type = "ATOM"
        
        # Get chains to extract based on proximity to reference
        chains_to_extract = get_chains_to_extract(
            chain_map, 
            reference_chain_ids,
            record_type,
            fixer_original, 
            params['cutoff_angstrom'], 
        )
        
        print(f"• Chains to extract NEAR {record_type} CHAINS: {chains_to_extract}")


    #=======================================================================
    # STEP 3: BINDING SITE EXTRACTION
    #=======================================================================
    
    params['verbose'] and add_separator("STEP 3: EXTRACTING CHAINS TO PDB")
    
    if chains_to_extract is not None:    
        extract_chains_to_pdb(file_paths['cleaned'], file_paths['binding_site'], chains_to_extract)
    else:
        # Default: If no chains could be identified, use the entire structure
        file_paths['binding_site'] = file_paths['cleaned']
        params['verbose'] and print("No ligand chains specified, using the entire structure")
    

    #=======================================================================
    # STEP 4: MISSING RESIDUE ANALYSIS
    #=======================================================================
    
    params['verbose'] and add_separator("STEP 4: COUNTING MISSING RESIDUES")
    
    # Analyze the original structure for missing residues
    original_pdbfixer = PDBFixer(filename=file_paths['input'])
  
    
    # Get missing residues by chain and chain ID mapping
    found_missing_residues, found_missing_chain_id_to_index = get_missing_residues_by_chain(
        original_pdbfixer, 
        chains_to_extract, 
        include_n_terminal_gaps=params['include_terminal_gaps'], 
        verbose=params['verbose']
    )

    # Clean structure again to exclude any remaining HETATM
    save_clean_structure(
        file_paths['binding_site'], 
        file_paths['cleaned'], 
        skip_hetatm=True, 
        verbose=params['verbose']
    )
    

    #=======================================================================
    # STEP 5-6: STRUCTURE COMPLETION
    #=======================================================================
    
    if len(found_missing_residues) > 0:
        params['verbose'] and add_separator("STEP 5: COMPLETE MISSING RESIDUES & ATOMS")
        
        # Load cleaned structure
        clean_fixer_object = load_structure_as_pdbfixer(file_paths['cleaned'])

        # Complete structure with missing residues and atoms
        completed_refined_fixer = complete_missing_structure(
            clean_fixer_object, 
            found_missing_residues=found_missing_residues,
            found_missing_chain_id_to_index=found_missing_chain_id_to_index,
            verbose=params['verbose']
        )
        
        params['verbose'] and add_separator("STEP 6: SAVE FINAL (COMPLETE) STRUCTURE")
        
        # Save the completed structure
        save_fixer_structure_to_pdb(
            completed_refined_fixer, 
            file_paths['temp_refined'], 
            verbose=params['verbose']
        )

        # PDBFixer renames chains - restore original chain IDs
        restore_original_chain_ids(
            file_paths['temp_refined'], 
            file_paths['refined'], 
            chains_to_extract, 
            verbose=params['verbose']
        )
    else:
        params['verbose'] and add_separator("STEP 5: SKIPPING STRUCTURE COMPLETION")
        
        # No missing residues, just restore chain IDs
        restore_original_chain_ids(
            file_paths['cleaned'], 
            file_paths['refined'], 
            chains_to_extract, 
            verbose=params['verbose']
        )
    

    #=======================================================================
    # STEP 7-11: STRUCTURE OPTIMIZATION AND CONVERSION
    #=======================================================================
   
   # Decide between fast track (skipping MolProbity) or full optimization
    if params['skip_molprobity']:
        # FAST TRACK: Skip MolProbity and directly use PDB2PQR
        params['verbose'] and add_separator("STEP 7: SKIPPING MOLPROBITY")
        params['verbose'] and add_separator("STEP 8: PROTONATION WITH PDB2PQR")
        
        # Protonate with PDB2PQR
        run_program(
            "PDB2PQR", 
            file_paths['refined'], 
            file_paths['protonated_pqr'], 
            verbose=params['verbose'], 
            pH_value=params['ph_value'], 
            config_file=params['config_file']
        )
        
        # Create PDBQT file for docking
        params['verbose'] and add_separator("STEP 9: CREATING PDBQT FILE FOR AutoDock Vina")  
        run_program(
            "MGLTools", 
            file_paths['protonated_pqr'], 
            file_paths['docking'], 
            verbose=params['verbose'], 
            pH_value=params['ph_value'], 
            config_file=params['config_file']
        )

        # params['verbose'] and add_separator("STEP 10: VALIDATING PDBQT FILE")  
        # validate_pdbqt_file(
        #     file_paths['docking'],
        #     file_paths['docking_final'],
        #     ['B'],
        #     verbose=params['verbose']
        # )

    else:
        # FULL OPTIMIZATION: Use MolProbity for hydrogen optimization
        
        #Step 1: Optimize structure with MolProbity
        params['verbose'] and add_separator("STEP 7: OPTIMIZING STRUCTURE FOR DOCKING")
        run_program(
            "MolProbity", 
            file_paths['refined'], 
            file_paths['optimized'], 
            verbose=params['verbose'], 
            pH_value=params['ph_value'], 
            config_file=params['config_file']
        )
        
        # Step 2: Process MolProbity output with OpenBabel
        params['verbose'] and add_separator("STEP 8: POSTPROCESSING MOLPROBITY FILE")   
        run_program(
            "OpenBabel", 
            file_paths['optimized'], 
            file_paths['edited'], 
            verbose=params['verbose'], 
            pH_value=params['ph_value'], 
            config_file=params['config_file']
        )

        fix_pdb_for_pdb2pqr(file_paths['edited'], file_paths['edited_no_conect'])

        # Step 3: Protonate with PDB2PQR
        params['verbose'] and add_separator("STEP 9: PROTONATION WITH PDB2PQR")
        run_program(
            "PDB2PQR", 
            file_paths['edited_no_conect'], 
            file_paths['protonated_pqr'], 
            verbose=params['verbose'], 
            pH_value=params['ph_value'], 
            config_file=params['config_file']
        )

        # Step 4: Create PDBQT file for AutoDock Vina
        params['verbose'] and add_separator("STEP 10: CREATING PDBQT FILE FOR AutoDock Vina")  
        run_program(
            "MGLTools", 
            file_paths['protonated_pqr'], 
            file_paths['docking'], 
            verbose=params['verbose'], 
            pH_value=params['ph_value'], 
            config_file=params['config_file']
        )

        # params['verbose'] and add_separator("STEP 11: VALIDATING PDBQT FILE")  
        # validate_pdbqt_file(
        #     file_paths['docking'],
        #     file_paths['docking'],
        #     ['B'],
        #     verbose=params['verbose']
        # )
    

    #=======================================================================
    # STEP 12: RESULTS SUMMARY
    #=======================================================================

    print("\nOutput files:")
    print(f"  {'File Type':<30} {'Path':<50}")
    print(f"  {'-'*30} {'-'*50}")
    print(f"  {'Binding site structure':<30} {os.path.relpath(file_paths['binding_site']):<50}")
    print(f"  {'Refined structure':<30} {os.path.relpath(file_paths['refined']):<50}")
    print(f"  {'Optimized structure':<30} {os.path.relpath(file_paths['edited']):<50}")
    print(f"  {'Final PDBQT for docking':<30} {os.path.relpath(file_paths['docking']):<50}")
    print("\n✅ Structure preparation completed successfully!")
    
    return file_paths


def main():
    """
    Execute the main workflow with command line arguments.
    
    This function parses command line arguments and executes the structure preparation workflow.
    """
    # Parse command line arguments
    args = parse_arguments()
    
    # Print verbose mode status
    args.verbose and print("\nRunning in verbose mode with detailed output...")
    
    #Execute the main workflow
    run_pipeline(args)


if __name__ == "__main__":
    main()