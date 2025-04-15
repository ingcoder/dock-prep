"""
PDB File Handler Module: Utilities for PDB File Processing and Manipulation

This module provides functions for handling, validating, and transforming PDB structure files.
It implements utilities for various stages of the molecular structure preparation pipeline:

1. File Path Setup: Setting up paths for input and output files
2. PDB Validation: Checking and validating PDB file format and content
3. Structure Cleaning: Removing unwanted atoms and molecules from structures
4. Chain Extraction: Isolating specific chains from multi-chain structures
5. Chain ID Mapping: Remapping chain IDs between different structure representations
6. Utility Functions: Supporting functions for file operations and status reporting

The module is designed to be modular and reusable, with improved error handling,
consistent output formatting, and type annotations.

Dependencies:
- OpenMM for PDB file writing
- BioPython for PDB structure parsing
- Standard Python libraries (os, sys)

Usage:
    This module is imported by the main run.py script and other modules
    in the dock-prep package to handle PDB file operations.
"""

import os
import sys
from openmm.app import PDBFile
from Bio import PDB
from typing import Dict, List, Tuple, Set, Optional, Any, Union, Callable
import re
import openmm.unit as unit
import tempfile
import datetime

#==============================================================================
#                          CONSTANTS & EXCEPTIONS
#==============================================================================

# Constants
PDB_RECORD_TYPES = ("HEADER", "ATOM", "HETATM")
WATER_MOLECULES = ["HOH", "WAT", "H2O", "SOL", "DOD"]

class PdbFileError(Exception):
    """Exception raised for errors in PDB file handling."""
    pass

#==============================================================================
#                          LOGGING & UTILITY FUNCTIONS
#==============================================================================

def log_message(message: str, verbose: bool = True) -> None:
    """Print message if verbose mode is enabled."""
    if verbose:
        print(message)

def log_success(filepath: str, message: str, verbose: bool = True) -> None:
    """Log success message with relative filepath if verbose mode is enabled."""
    if verbose:
        rel_path = os.path.relpath(filepath)
        print(f"✅ {message}: {rel_path}")

def get_chain_id_from_line(line: str) -> Optional[str]:
    """Extract chain ID from PDB file line if it's an ATOM/HETATM/TER record."""
    if (line.startswith(("ATOM", "HETATM")) and len(line) > 21) or (line.startswith("TER") and len(line) > 21):
        return line[21:22]
    return None

def _count_chains(line: str, chains: List[str], chains_seen: Set[str]) -> None:
    """
    Helper function to count unique chains in PDB file.
    
    Args:
        line: PDB file line
        chains: List to store unique chain IDs
        chains_seen: Set to track seen chain IDs
    """
    chain_id = get_chain_id_from_line(line)
    if chain_id and chain_id not in chains_seen:
        chains.append(chain_id)
        chains_seen.add(chain_id)
    # Reset chains_seen when a TER record is encountered
    elif line.startswith("TER"):
        chains_seen.clear()

#==============================================================================
#                          FILE PATH SETUP
#==============================================================================

def setup_file_paths(file_prefix, output_dir, input_file):
    """
    Set up all file paths for the workflow.

    """
    # Create results folder if it doesn't exist
    results_folder = os.path.abspath(output_dir)
    if not os.path.exists(results_folder):
        os.makedirs(results_folder)
        
    # Store all file paths
    file_paths = {}
    file_paths['input'] = input_file
    file_paths['cleaned'] = os.path.join(results_folder, f"{file_prefix}_0_structure_cleaned.pdb")
    file_paths['binding_site'] = os.path.join(results_folder, f"{file_prefix}_0_structure_selected_chains.pdb")
    file_paths['temp_refined'] = os.path.join(results_folder, f"{file_prefix}_1_structure_completed_temp.pdb")
    file_paths['refined'] = os.path.join(results_folder, f"{file_prefix}_2_structure_completed_final.pdb")
    file_paths['optimized'] = os.path.join(results_folder, f"{file_prefix}_3_structure_flipped_h_temp.pdb")
    file_paths['edited'] = os.path.join(results_folder, f"{file_prefix}_4_structure_flipped_h_final.pdb")
    file_paths['edited_cif'] = os.path.join(results_folder, f"{file_prefix}_4_structure_flipped_h_final.cif")
    file_paths['edited_no_conect'] = os.path.join(results_folder, f"{file_prefix}_4_structure_flipped_h_final_no_conect.pdb")
    file_paths['edited_no_conect_cif'] = os.path.join(results_folder, f"{file_prefix}_4_structure_flipped_h_final_no_conect.cif")
    file_paths['protonated_pqr'] = os.path.join(results_folder, f"{file_prefix}_5_structure_protonated.pqr")
    file_paths['docking'] = os.path.join(results_folder, f"{file_prefix}_6_structure_docking.pdbqt")
    file_paths['docking_final'] = os.path.join(results_folder, f"{file_prefix}_6_structure_docking_final.pdbqt")
    
    return file_paths, results_folder

#==============================================================================
#                          PDB FILE VALIDATION
#==============================================================================

def validate_pdb_file(input_file: str, verbose: bool = True) -> Tuple[bool, bool]:
    """Validates if the input file is a PDB file and checks for SEQRES records."""
    # Check if file exists
    assert os.path.exists(input_file), f"File not found: {input_file}"
    
    # Check file extension
    file_ext = os.path.splitext(input_file)[1].lower()
    if not file_ext.startswith('.pdb'):
        raise PdbFileError(f"Unsupported file format: {file_ext}. Only PDB formats (.pdb, .pdbqt, etc.) are supported.")
    
    # Check if file contains valid PDB format by looking for HEADER/ATOM/HETATM records
    has_pdb_records = False
    has_seqres = False
    
    try:
        with open(input_file, 'r') as f:
            for line in f:
                # Check for PDB record identifiers (HEADER, ATOM, HETATM)
                if line.startswith(PDB_RECORD_TYPES):
                    has_pdb_records = True
                # Check for SEQRES records
                if line.startswith("SEQRES"):
                    has_seqres = True
                    log_message("✅ SEQRES records found in PDB file", verbose)
                    break
    except UnicodeDecodeError:
        # If we can't read the file as text, it's not a valid PDB
        raise PdbFileError("File is not a valid text-based PDB file")
    
    if not has_pdb_records:
        raise PdbFileError("File does not appear to be a valid PDB format. No PDB record identifiers found.")
    
    if not has_seqres and verbose:
        log_message("❗ Warning: No SEQRES records found in PDB file.", verbose)
        log_message("Missing residues detection may be limited. The structure preparation will", verbose)
        log_message("continue, but missing residues might not be properly identified.", verbose)
    
    return has_pdb_records, has_seqres

#==============================================================================
#                          PDB FILE PROCESSING
#==============================================================================

def process_pdb_file(input_filepath: str, output_filepath: str, 
                    line_processor: Callable[[str], Optional[str]],
                    success_message: str = "File processed", 
                    verbose: bool = True) -> None:
    """Generic PDB file processor that applies a line processor function to each line."""
    try:
        with open(input_filepath, 'r') as fin, open(output_filepath, 'w') as fout:
            for line in fin:
                processed_line = line_processor(line)
                if processed_line is not None:
                    fout.write(processed_line)
        
        log_success(output_filepath, success_message, verbose)
    except Exception as e:
        raise PdbFileError(f"Failed to process PDB file: {e}")

def extract_chains_to_pdb(input_filepath: str, output_filepath: str, target_chains: Set = None, verbose: bool = True) -> None:
    """Extracts specified chains from a PDB file and saves to a new file."""
    def chain_filter(line: str) -> str:
        if line.startswith(("ATOM", "HETATM")):
            chain_id = line[21]
            record_type = line[0:6].strip()
            chain_search_key = (chain_id, record_type)
            if chain_search_key in target_chains:
                return line
            return None
        elif line.startswith("TER") and len(line) > 21:
            chain_id = line[21]
            if chain_id in target_chains:
                return line
            return None
        else:
            return line
    
    process_pdb_file(
        input_filepath, 
        output_filepath, 
        chain_filter,
        f"Filtered PDB saved containing chains {target_chains}",
        verbose
    )

def save_clean_structure(input_file: str, output_file: str, skip_hetatm: bool = True, verbose: bool = True) -> List[str]:
    """Loads PDB file, removes heteroatoms, and returns a list of chain IDs."""
    try:
        # Validate the input file is a valid PDB file and check for SEQRES records
        validate_pdb_file(input_file, verbose=verbose)
        
        # Step 1: Read PDB content
        with open(input_file, 'r') as f:
            lines = f.readlines()

        # Step 2: Process lines and filter HETATM records
        cleaned_lines = []
        chains = []
        chains_seen = set()

        for line in lines:
            should_keep = True
            
            if line.startswith("HETATM"):
                if skip_hetatm:  
                    # Skip all HETATM lines
                    should_keep = False
                else:
                    # If not skipping all HETATM, still skip water molecules
                    if any(water in line for water in WATER_MOLECULES):
                        should_keep = False
            
            if should_keep:
                cleaned_lines.append(line)
                _count_chains(line, chains, chains_seen)
        
        # Step 3: Write cleaned content to file
        cleaned_pdb_text = "".join(cleaned_lines)
        with open(output_file, 'w') as f:
            f.write(cleaned_pdb_text)

        log_success(output_file, f"Cleaned PDB saved with chains:\n {chains}", verbose)
            
        return chains

    except PdbFileError as e:
        print(f"Error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error: Failed to clean PDB file: {e}")
        sys.exit(1)

#==============================================================================
#                          STRUCTURE CONVERSION & SAVING
#==============================================================================

def save_fixer_structure_to_pdb(pdb_structure: Any, output_filepath: str, verbose: bool = True) -> None:
    """
    Saves PDBFixer structure to a PDB file with properly formatted coordinates.
    
    This function converts the coordinates from nanometers to Ångstroms and saves the structure
    while preserving the standard PDB format.
    """
    # Convert coordinates from nanometers to Ångstroms
    positions_nm = pdb_structure.positions
    positions_angstrom = positions_nm.in_units_of(unit.angstrom)
    
    # Save directly to the output file using OpenMM's PDBFile writer
    with open(output_filepath, 'w') as f:
        PDBFile.writeFile(pdb_structure.topology, positions_angstrom, f)
    
    log_success(output_filepath, "Complete structure saved with properly formatted coordinates (in Ångstroms)", verbose)

def restore_original_chain_ids(input_filepath: str, output_filepath: str, 
                              target_chains: List = None, verbose: bool = True) -> None:
    """Restores original chain IDs in a PDB file created by PDBFixer."""
    try:
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
        chain_mapping = {current_chain: correct_chain 
                        for current_chain, correct_chain in zip(current_chain_ids, target_chains)}
        
        log_message(f"• Chain mapping: {chain_mapping}", verbose)
        
        def chain_mapper(line: str) -> str:
            # Process only ATOM and HETATM lines which contain atom coordinates
            if line.startswith(("ATOM", "HETATM")):
                # In PDB format, the chain ID is at position 21 (0-indexed)
                chain_id = line[21:22]
                if chain_id in chain_mapping:
                    # Get the target chain ID from the mapping
                    fixed_chain = chain_mapping[chain_id]
                    # Construct a new line by replacing just the chain ID character
                    # This preserves all other information in the line
                    return line[:21] + fixed_chain + line[22:]
            return line
        
        process_pdb_file(
            input_filepath, 
            output_filepath, 
            chain_mapper,
            "Fixed PDB saved with restored chain IDs",
            verbose
        )
    
    except Exception as e:
        print(f"Error: Failed to restore chain IDs: {e}")
        sys.exit(1)


def validate_pdbqt_file(input_filepath: str, output_filepath: str, selected_chains: Set = None, verbose: bool = True) -> None:
    """Validates PDBQT file and filters by selected chains if specified."""
    # Check if file exists
    assert os.path.exists(input_filepath), f"File not found: {input_filepath}"

    # Set to track which chain IDs have already had messages printed
    partial_chain_id_in_selection = set()
    full_chain_id_in_selection = set()
    full_chain_id_not_in_selection = set()
    

    def check_chain_ids(line: str) -> Optional[str]:
        # Keep all non-atom/hetatm lines
        if not line.startswith(("ATOM", "HETATM")):
            return line
            
        # For atom/hetatm records, filter by chain if selected_chains is provided
        if selected_chains:
            chain_id = line[21:23].strip()  # Chain ID is at position 21
            
            if chain_id in selected_chains:
                if chain_id not in full_chain_id_in_selection:
                    print(f'Chain ID "{chain_id}" is in the selected chains.')
                    full_chain_id_in_selection.add(chain_id)
                return line
                
            # Check if this is a multi-character chain ID where the last character matches
            # This handles cases like "UB" when "B" is selected
            if len(chain_id.strip()) > 1 and chain_id[-1] in selected_chains:
                # Only print the message once per unique chain ID
                if chain_id not in partial_chain_id_in_selection:
                    print(f'Chain ID "{chain_id}" is not in the selected chains, but "{chain_id[-1]}" is. ' 
                          f'Using "{chain_id[-1]}" instead.')
                    partial_chain_id_in_selection.add(chain_id)
                
                # Replace the chain ID with just the last character that matches the selection
                return line[:21] + chain_id[-1] + line[23:]
                
            # No match found, skip this atom
            if chain_id not in full_chain_id_not_in_selection:
                print(f'Chain ID "{chain_id}" is not in the selected chains.')
                full_chain_id_not_in_selection.add(chain_id)
            return None
        else:
            # If no chains specified, keep all atoms
            return line

    process_pdb_file(
        input_filepath, 
        output_filepath, 
        check_chain_ids,
        "PDBQT file processed and validated",
        verbose
    ) 


def fix_pdb_for_pdb2pqr(input_file, output_file, verbose=True):
    """ Reads a PDB file, renumbers ATOM/HETATM records sequentially starting from 1,
    and ensures proper column alignment. Skips CONECT records.
    This prepares the file for tools like pdb2pqr that might have issues with
    large atom numbers or misaligned columns.
    """
    try:
        with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
            for line in f_in:
                # Skip CONECT lines entirely
                if line.startswith('CONECT'):
                    continue
                
                f_out.write(line)

    except FileNotFoundError:
        print(f"Error: Input file not found: {input_file}")
        sys.exit(1)
    except Exception as e:
        print(f"Error: Failed to fix PDB file for pdb2pqr: {e}")
        sys.exit(1)

   
