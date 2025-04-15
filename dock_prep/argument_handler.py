import argparse
import os

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
    parser.add_argument("--file_prefix", type=str, default="prefix", help="Example PDB identifier (e.g., 1JY7)")
    parser.add_argument("--include_terminal_gaps", action="store_true", help="Include N-terminal gaps for missing residues analysis")
    
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

def extract_params_from_args(args):
    """
    Extract and process parameters from command line arguments.
    
    Args:
        args: Command line arguments parsed by argparse
        
    Returns:
        dict: Dictionary containing all extracted parameters
    """
    params = {}
    
    # Extract basic parameters
    params['input_file'] = args.input_file
    params['verbose'] = args.verbose
    params['file_prefix'] = args.file_prefix
    params['cutoff_angstrom'] = args.cutoff
    params['results_folder'] = args.output_dir
    params['skip_molprobity'] = args.skip_molprobity
    params['ph_value'] = args.ph
    params['include_terminal_gaps'] = args.include_terminal_gaps
    
    # Parse chain IDs
    params['include_chain_ids'] = args.include_chains.split(',') if args.include_chains else None
    params['reference_atom_chain_ids'] = args.reference_atom_chains.split(',') if args.reference_atom_chains else None
    params['reference_hetatm_chain_ids'] = args.reference_hetatm_chains.split(',') if args.reference_hetatm_chains else None
    
    # Setup config file path
    script_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    params['config_file'] = os.path.join(script_dir, 'scripts', 'config_env.json')

    # Determine whether to skip HETATM records
    params['skip_hetatm'] = True if params['reference_hetatm_chain_ids'] is None else False   



    # Print summary of output files
    print("\nSettings:")
    print(f"  {'Parameter':<30} {'Value':<50}")
    print(f"  {'-'*30} {'-'*50}")
    print(f"  {'INPUT_FILE':<30} {params['input_file']}")
    print(f"  {'RESULTS_FOLDER':<30} {params['results_folder']}")
    print(f"  {'CONFIG_FILE':<30} {params['config_file']}")
    print(f"  {'VERBOSE':<30} {params['verbose']}")
    print(f"  {'FILE_PREFIX':<30} {params['file_prefix']}")
    print(f"  {'INCLUDE_CHAIN_IDS':<30} {params['include_chain_ids']}")
    print(f"  {'REFERENCE_ATOM_CHAIN_IDS':<30} {params['reference_atom_chain_ids']}")
    print(f"  {'REFERENCE_HETATM_CHAIN_IDS':<30} {params['reference_hetatm_chain_ids']}")
    print(f"  {'CUTOFF_ANGSTROM':<30} {params['cutoff_angstrom']}")
    print(f"  {'pH_VALUE':<30} {params['ph_value']}")
    print(f"  {'SKIP_HETATM':<30} {params['skip_hetatm']}")
    return params 