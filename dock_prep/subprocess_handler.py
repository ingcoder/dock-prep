"""
Subprocess handler for molecular structure file conversions.

This module provides utilities to run external programs like MGLTools, OpenBabel,
MolProbity, and PDB2PQR for converting between different molecular structure formats,
with a focus on preparing PDBQT files for molecular docking. It constructs and runs
shell commands to interface with these external tools.
"""

import subprocess
import os
import traceback
import json

def _check_if_file_exists(input_path):
    """
    Checks if a file exists and returns the absolute path.
    """
    if not os.path.exists(input_path):
        print(f"⚠️ WARNING: Input file not found at: {input_path}")
        # List directory contents to see what's there
        directory = os.path.dirname(input_path) or '.'
        print(f"Contents of directory {directory}:")
        for file in os.listdir(directory):
            print(f"  - {file}")
        raise FileNotFoundError(f"Input file not found: {input_path}")
    else:
        print(f"✅ Input file exists: {input_path}")
        print(f"File size: {os.path.getsize(input_path)} bytes")
        return True

def get_conda_env():
    """
    Returns the name of the current conda environment.
    """
    conda_env = os.environ.get('CONDA_DEFAULT_ENV', None)
    if conda_env:
        return conda_env
    else:
        # Alternative method if CONDA_DEFAULT_ENV is not set
        try:
            result = subprocess.run('conda info --envs | grep "*"', 
                                   shell=True, 
                                   capture_output=True, 
                                   text=True)
            if result.returncode == 0 and result.stdout:
                # Parse the output to get the environment name
                return result.stdout.strip().split()[0]
        except Exception:
            pass
        return "No conda environment detected"
    
def get_env_vars(config_file):
    """
    Returns the environment variables for a given tool.
    """
    if not os.path.exists(config_file):
        raise FileNotFoundError(f"Environment file not found at: {config_file}")
    with open(config_file, 'r') as f:
        env_vars = json.load(f)
    return env_vars

def _run_subprocess_command(tool_name, command, abs_input_path, abs_output_path, verbose=True, config_file=None):
    """
    Runs a subprocess command and returns the output.
    """
    # Run the command using subprocess.run
    print(f"Running command: {tool_name}")
 
    # Get the current working directory to restore it later
    original_cwd = os.getcwd()
    print(f"Original working directory: {original_cwd}")
    
    # Print the conda environment
    conda_env = get_conda_env()
    print(f"Current conda environment: {conda_env}")
    
    try:
        # MGLTools script is looking for the file in the current working directory,
        # not in the absolute path you're providing, so we need to change the working directory
        if tool_name == "MGLTools":
            input_dir = os.path.dirname(abs_input_path)
            input_filename = os.path.basename(abs_input_path)
            output_filename = os.path.basename(abs_output_path)
            print(f"input_dir: {input_dir}")
            
            # Change working directory
            if input_dir:
                print(f"Changing directory to: {input_dir}")
                os.chdir(input_dir)
            
            # Update command to use relative file paths - use named parameters to avoid confusion
            command = construct_shell_command(
                tool_name=tool_name, 
                abs_input_path=input_filename, 
                abs_output_path=output_filename, 
                pH_value=7.4,  # Default pH value - consider passing this from the calling function if needed
                config_file=config_file
            )
            # verbose and print(f"Updated command: {command}")
        
        # Run the command
        result = subprocess.run(command, shell=True, capture_output=True, text=True, timeout=300, check=True)
        
        # # Print command output if verbose
        # if verbose and result.stdout:
        #     print("Command output:")
        #     print(result.stdout)
        
        # Check if the output file was created
        if os.path.exists(abs_output_path):
            print(f"Successfully created output file: {abs_output_path}")
            print(f"File size: {os.path.getsize(abs_output_path)} bytes")
            return True
        else:
            print("Error: Output file was not created despite successful command execution")
            return False
    except subprocess.CalledProcessError as e:
        print(f"Command execution failed with return code {e.returncode}")
        if e.stdout:
            print("Standard output:")
            print(e.stdout)
        if e.stderr:
            print("Standard error:")
            print(e.stderr)
        raise
    except Exception as e:
        print(f"Unexpected error: {str(e)}")
        traceback.print_exc()
        raise e
    finally:
        # Always restore the original working directory
        if os.getcwd() != original_cwd:
            os.chdir(original_cwd)
            # verbose and print(f"Restored directory to: {original_cwd}")

def construct_shell_command(tool_name, abs_input_path, abs_output_path, pH_value, config_file=None):
    """
    Constructs a shell command for a given tool and command.
    """
    # Extract file extensions
    input_ext = os.path.splitext(abs_input_path)[1].lower().lstrip('.')
    output_ext = os.path.splitext(abs_output_path)[1].lower().lstrip('.')
    
    if tool_name == "MGLTools":
        # Load MGLTools environment variables
        env_vars = get_env_vars(config_file)
        pythonsh = os.path.join(env_vars['MGL_BIN'], 'pythonsh')
        prepare_receptor_script = os.path.join(env_vars['MGL_PACKAGES'], 'AutoDockTools', 'Utilities24', 'prepare_receptor4.py')
        
        _check_if_file_exists(pythonsh)
        _check_if_file_exists(prepare_receptor_script)
   
        # MGLTools has issues with absolute paths, so use the filename directly
        input_filename = os.path.basename(abs_input_path) if os.path.isabs(abs_input_path) else abs_input_path
        output_filename = os.path.basename(abs_output_path) if os.path.isabs(abs_output_path) else abs_output_path
        
        # The "checkhydrogens" option tells the script to "add hydrogens only if there are none already" in the structure.
        # -C flag tells the script to preserve input charges. It will not add new Gasteiger charges
        # -U nphs_lps flag. - U defines type of cleanup. nphs_lps remove non-polar hydrogen atoms (C-H hydrogens) and merge their charges with their attached carbon atoms
        # It's worth noting that the default cleanup setting is more extensive (nphs_lps_waters_nonstdres) Removes water and non standard residues.

        if input_ext == "pqr": 
            return f"conda run -n {env_vars['MGL_ENV_NAME']} --no-capture-output bash -c 'export PYTHONPATH={env_vars['MGL_PACKAGES']}:$PYTHONPATH && {pythonsh} {prepare_receptor_script} -r {input_filename} -o {output_filename} -A checkhydrogens -C -U nphs_lps'"
        else:
            return f"conda run -n {env_vars['MGL_ENV_NAME']} --no-capture-output bash -c 'export PYTHONPATH={env_vars['MGL_PACKAGES']}:$PYTHONPATH && {pythonsh} {prepare_receptor_script} -r {input_filename} -o {output_filename} -A checkhydrogens'"
        # Command for prepare_receptor4.py using conda run (simpler approach)
        # return f"conda run -n {env_vars['MGL_ENV_NAME']} {pythonsh} {prepare_receptor_script} -r {input_filename} -o {output_filename} -A checkhydrogens"
    elif tool_name == "OpenBabel":
        # Map file extensions to OpenBabel format codes
        format_map = {'pdb': 'pdb','pqr': 'pqr','mol': 'mol','mol2': 'mol2','sdf': 'sdf','xyz': 'xyz','pdbqt': 'pdbqt','mmcif': 'cif','cif': 'cif'}
        input_format = format_map.get(input_ext)
        output_format = format_map.get(output_ext)
        
        # Return the correct OpenBabel command with proper format specifiers
        return f"obabel -i{input_format} {abs_input_path} -o{output_format} -O {abs_output_path}"
    elif tool_name == "MolProbity":
        # Define paths to MolProbity executables
        env_vars = get_env_vars(config_file)
        reduce = os.path.join(env_vars['MOLPROBITY_BIN'],'reduce')
        probe = os.path.join(env_vars['MOLPROBITY_BIN'],'probe')
        _check_if_file_exists(reduce)
        _check_if_file_exists(probe)

        return f"{reduce} -FLIP {abs_input_path} > {abs_output_path}"
    elif tool_name == "PDB2PQR":
        return f"pdb2pqr30 --ff AMBER --keep-chain --titration-state-method propka --with-ph {pH_value} {abs_input_path} {abs_output_path}"
   
def run_program(tool_name, input_path, output_path, verbose=True, pH_value=7.4, config_file=None):
    """
    Runs a program with the given tool name and command.
    """
    # Convert to absolute paths
    abs_input_path = os.path.abspath(input_path)
    abs_output_path = os.path.abspath(output_path)

    _check_if_file_exists(abs_input_path)
    
    # Use named parameters to avoid confusion with positional parameters
    command = construct_shell_command(tool_name=tool_name, 
                               abs_input_path=abs_input_path, 
                               abs_output_path=abs_output_path, 
                               pH_value=pH_value,
                               config_file=config_file)
    
    # Check if a valid command was returned
    if not command:
        raise ValueError(f"Failed to construct a valid command for {tool_name}")
        
    _run_subprocess_command(tool_name, command, abs_input_path, abs_output_path, verbose, config_file)