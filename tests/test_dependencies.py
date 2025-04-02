import unittest
import importlib
import subprocess
import os
import shutil
import json

class TestDependencies(unittest.TestCase):
    """Tests to verify all required dependencies are properly installed."""
    
    def setUp(self):
        """Set up test environment."""
        # Try to load the config file if it exists
        self.config_file = None
        project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        possible_config = os.path.join(project_root, "dock_prep", "config.json")
        if os.path.exists(possible_config):
            self.config_file = possible_config
    
    def test_python_dependencies(self):
        """Test that all required Python packages are installed."""
        required_packages = [
            'numpy',
            'pdbfixer',
            'openmm',
            'Bio' #installed with conda install biopython
        ]
        
        missing_packages = []
        for package in required_packages:
            try:
                importlib.import_module(package)
                print(f"✅ Package {package} is installed")
            except ImportError:
                missing_packages.append(package)
                print(f"❌ Package {package} is NOT installed")
        
        self.assertEqual(len(missing_packages), 0, 
                         f"The following required packages are missing: {', '.join(missing_packages)}")
    
    def test_conda_environment(self):
        """Test that we're in a conda environment."""
        from dock_prep.subprocess_handler import get_conda_env
        
        conda_env = get_conda_env()
        print(f"Current conda environment: {conda_env}")
        self.assertNotEqual(conda_env, "No conda environment detected", 
                           "Not running in a conda environment. This package requires conda.")
    
    def test_external_tools_existence(self):
        """Test that external tools are installed and available on the PATH."""
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
        
        self.assertEqual(len(missing_tools), 0, 
                         f"The following required tools are missing from PATH: {', '.join(missing_tools)}")
    
    def test_mgltools_installation(self):
        """Test that MGLTools is properly configured."""
        if not self.config_file:
            self.skipTest("Config file not found, skipping MGLTools test")
        
        try:
            with open(self.config_file, 'r') as f:
                config = json.load(f)
            
            # Check if MGLTools environment variables are set
            self.assertIn('MGL_BIN', config, "MGL_BIN not found in config file")
            self.assertIn('MGL_PACKAGES', config, "MGL_PACKAGES not found in config file")
            self.assertIn('MGL_ENV_NAME', config, "MGL_ENV_NAME not found in config file")
            
            # Check if pythonsh exists
            pythonsh = os.path.join(config['MGL_BIN'], 'pythonsh')
            self.assertTrue(os.path.exists(pythonsh), f"MGLTools pythonsh not found at: {pythonsh}")
            
            # Check if prepare_receptor4.py exists
            prepare_receptor = os.path.join(config['MGL_PACKAGES'], 'AutoDockTools', 
                                           'Utilities24', 'prepare_receptor4.py')
            self.assertTrue(os.path.exists(prepare_receptor), 
                           f"prepare_receptor4.py not found at: {prepare_receptor}")
            
            print(f"✅ MGLTools is properly configured")
        except (FileNotFoundError, json.JSONDecodeError, KeyError) as e:
            self.fail(f"MGLTools configuration check failed: {str(e)}")
    
    def test_molprobity_installation(self):
        """Test that MolProbity is properly configured."""
        if not self.config_file:
            self.skipTest("Config file not found, skipping MolProbity test")
        
        try:
            with open(self.config_file, 'r') as f:
                config = json.load(f)
            
            # Check if MolProbity environment variable is set
            self.assertIn('MOLPROBITY_BIN', config, "MOLPROBITY_BIN not found in config file")
            
            # Check if reduce and probe exist
            reduce = os.path.join(config['MOLPROBITY_BIN'], 'reduce')
            probe = os.path.join(config['MOLPROBITY_BIN'], 'probe')
            
            self.assertTrue(os.path.exists(reduce), f"MolProbity reduce not found at: {reduce}")
            self.assertTrue(os.path.exists(probe), f"MolProbity probe not found at: {probe}")
            
            print(f"✅ MolProbity is properly configured")
        except (FileNotFoundError, json.JSONDecodeError, KeyError) as e:
            self.fail(f"MolProbity configuration check failed: {str(e)}")


if __name__ == '__main__':
    unittest.main() 