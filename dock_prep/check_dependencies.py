"""
Utility module for checking if all required dependencies are correctly installed.
"""

import importlib
import os
import shutil
import json
import sys

def check_python_packages():
    """
    Check if all required Python packages are installed.
    
    Returns:
        tuple: (all_installed, missing_packages)
    """
    required_packages = [
        'numpy',
        'pdbfixer',
        'openmm',
        'Bio',
    ]
    
    missing_packages = []
    for package in required_packages:
        try:
            importlib.import_module(package)
            print(f"✅ Package {package} is installed")
        except ImportError:
            missing_packages.append(package)
            print(f"❌ Package {package} is NOT installed")
    
    return (len(missing_packages) == 0, missing_packages)

def check_conda_environment():
    """
    Check if running in a conda environment.
    
    Returns:
        bool: True if running in a conda environment, False otherwise
    """
    from dock_prep.subprocess_handler import get_conda_env
    
    conda_env = get_conda_env()
    print(f"Current conda environment: {conda_env}")
    
    if conda_env == "No conda environment detected":
        print("❌ Not running in a conda environment. This package requires conda.")
        return False
    return True

def check_external_tools():
    """
    Check if required external tools are available on PATH.
    
    Returns:
        tuple: (all_installed, missing_tools)
    """
    external_tools = {
        "OpenBabel": "obabel",
        "PDB2PQR": "pdb2pqr30"
    }
    
    missing_tools = []
    for tool_name, command in external_tools.items():
        if shutil.which(command) is None:
            missing_tools.append(tool_name)
            print(f"❌ Tool {tool_name} ({command}) is NOT on PATH")
        else:
            print(f"✅ Tool {tool_name} ({command}) is available")
    
    return (len(missing_tools) == 0, missing_tools)

def check_config_tools():
    """
    Check if MGLTools and MolProbity are properly configured.
    
    Returns:
        tuple: (has_config, mgltools_ok, molprobity_ok, config_errors)
    """
    # Try to load the config file if it exists
    config_errors = []
    config_file = None
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    possible_config = os.path.join(project_root, "scripts", "config_env.json")
    if os.path.exists(possible_config):
        config_file = possible_config
    else:
        print("❌ Config file not found at", possible_config)
        return (False, False, False, ["Config file not found"])
    
    try:
        with open(config_file, 'r') as f:
            config = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError) as e:
        print(f"❌ Error reading config file: {str(e)}")
        return (False, False, False, [f"Error reading config file: {str(e)}"])
    
    # Check MGLTools
    mgltools_ok = True
    # Check if MGLTools environment variables are set
    for var in ['MGL_BIN', 'MGL_PACKAGES', 'MGL_ENV_NAME']:
        if var not in config:
            print(f"❌ {var} not found in config file")
            config_errors.append(f"{var} missing from config")
            mgltools_ok = False
    
    if mgltools_ok:
        # Check if pythonsh exists
        pythonsh = os.path.join(config['MGL_BIN'], 'pythonsh')
        if not os.path.exists(pythonsh):
            print(f"❌ MGLTools pythonsh not found at: {pythonsh}")
            config_errors.append(f"MGLTools pythonsh not found at: {pythonsh}")
            mgltools_ok = False
        
        # Check if prepare_receptor4.py exists
        prepare_receptor = os.path.join(config['MGL_PACKAGES'], 'AutoDockTools', 
                                       'Utilities24', 'prepare_receptor4.py')
        if not os.path.exists(prepare_receptor):
            print(f"❌ prepare_receptor4.py not found at: {prepare_receptor}")
            config_errors.append(f"prepare_receptor4.py not found at: {prepare_receptor}")
            mgltools_ok = False
        
        if mgltools_ok:
            print(f"✅ MGLTools is properly configured")
    
    # Check MolProbity
    molprobity_ok = True
    # Check if MolProbity environment variable is set
    if 'MOLPROBITY_BIN' not in config:
        print(f"❌ MOLPROBITY_BIN not found in config file")
        config_errors.append("MOLPROBITY_BIN missing from config")
        molprobity_ok = False
    
    if molprobity_ok:
        # Check if reduce and probe exist
        reduce = os.path.join(config['MOLPROBITY_BIN'], 'reduce')
        probe = os.path.join(config['MOLPROBITY_BIN'], 'probe')
        
        if not os.path.exists(reduce):
            print(f"❌ MolProbity reduce not found at: {reduce}")
            config_errors.append(f"MolProbity reduce not found at: {reduce}")
            molprobity_ok = False
        
        if not os.path.exists(probe):
            print(f"❌ MolProbity probe not found at: {probe}")
            config_errors.append(f"MolProbity probe not found at: {probe}")
            molprobity_ok = False
        
        if molprobity_ok:
            print(f"✅ MolProbity is properly configured")
    
    return (True, mgltools_ok, molprobity_ok, config_errors)

def check_all_dependencies():
    """
    Check all dependencies and return a detailed report.
    
    Returns:
        bool: True if all dependencies are properly installed, False otherwise
    """
    print("="*50)
    print("DEPENDENCY CHECK REPORT")
    print("="*50)
    
    all_ok = True
    
    # Check Python packages
    print("\nChecking Python packages...")
    packages_ok, missing_packages = check_python_packages()
    if not packages_ok:
        all_ok = False
    
    # Check Conda environment
    print("\nChecking Conda environment...")
    conda_ok = check_conda_environment()
    if not conda_ok:
        all_ok = False
    
    # Check external tools
    print("\nChecking external tools...")
    tools_ok, missing_tools = check_external_tools()
    if not tools_ok:
        all_ok = False
    
    # Check config-based tools
    print("\nChecking configuration-based tools...")
    has_config, mgltools_ok, molprobity_ok, config_errors = check_config_tools()
    if not has_config or not mgltools_ok or not molprobity_ok:
        all_ok = False
    
    # Print summary
    print("\n" + "="*50)
    print("SUMMARY")
    print("="*50)
    
    if all_ok:
        print("✅ All dependencies are properly installed!")
    else:
        print("❌ Some dependencies are missing or improperly configured:")
        
        if not packages_ok:
            print(f"  - Missing Python packages: {', '.join(missing_packages)}")
        
        if not conda_ok:
            print("  - Not running in a conda environment")
        
        if not tools_ok:
            print(f"  - Missing external tools: {', '.join(missing_tools)}")
        
        if not has_config:
            print("  - Configuration file is missing or invalid")
        elif not mgltools_ok:
            print("  - MGLTools is not properly configured")
        elif not molprobity_ok:
            print("  - MolProbity is not properly configured")
        
        if config_errors:
            print("\nConfiguration errors:")
            for error in config_errors:
                print(f"  - {error}")
    
    return all_ok

if __name__ == "__main__":
    # When run as a script, check all dependencies and exit with appropriate code
    success = check_all_dependencies()
    sys.exit(0 if success else 1) 