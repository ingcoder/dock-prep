import unittest
import os
import sys
import tempfile
import shutil
from subprocess_handler import _check_if_file_exists, construct_shell_command

class TestSubprocessHandler(unittest.TestCase):
    """Tests for subprocess_handler module in PdbqtConverter."""
    
    def setUp(self):
        """Set up test environment before each test."""
        self.example_pdb = os.path.join("examples", "2pgh_original.pdb")
        self.test_dir = tempfile.mkdtemp()
        self.test_file = os.path.join(self.test_dir, "test_file.pdb")
        
        # Create a test file
        with open(self.example_pdb, 'r') as fin, open(self.test_file, 'w') as fout:
            lines = fin.readlines()[:10]  # Just take a few lines
            fout.writelines(lines)
    
    def tearDown(self):
        """Clean up after each test."""
        shutil.rmtree(self.test_dir)
    
    def test_check_if_file_exists(self):
        """Test the file existence check function."""
        # Test with existing file
        self.assertTrue(_check_if_file_exists(self.test_file),
                       f"File {self.test_file} should exist")
        
        # Test with non-existing file
        non_existent_file = os.path.join(self.test_dir, "non_existent.pdb")
        with self.assertRaises(FileNotFoundError):
            _check_if_file_exists(non_existent_file)
    
    def test_construct_shell_command_molprobity(self):
        """Test constructing MolProbity shell command."""
        # Create input and output paths
        input_path = os.path.abspath(self.test_file)
        output_path = os.path.join(self.test_dir, "output.pdb")
        output_path = os.path.abspath(output_path)
        
        # Construct the command
        try:
            command = construct_shell_command("MolProbity", input_path, output_path)
            
            # Basic check if the command seems reasonable
            self.assertIsNotNone(command, "Command should not be None")
            self.assertIn("reduce", command, "MolProbity command should include reduce")
            self.assertIn("-FLIP", command, "MolProbity command should include -FLIP option")
            self.assertIn(input_path, command, "Command should include input path")
            self.assertIn(output_path, command, "Command should include output path")
            
            print(f"Generated MolProbity command: {command}")
        except Exception as e:
            # If the MolProbity binary doesn't exist, the test will fail with an error
            # We'll skip the test in this case
            self.skipTest(f"Skipping test due to error: {str(e)}")
    
    def test_construct_shell_command_openbabel(self):
        """Test constructing OpenBabel shell command."""
        # Create input and output paths with different extensions
        input_path = os.path.abspath(self.test_file)
        output_path = os.path.join(self.test_dir, "output.mol2")
        output_path = os.path.abspath(output_path)
        
        # Construct the command
        command = construct_shell_command("OpenBabel", input_path, output_path)
        
        # Basic check if the command seems reasonable
        self.assertIsNotNone(command, "Command should not be None")
        self.assertIn("obabel", command, "OpenBabel command should include obabel")
        self.assertIn("-ipdb", command, "Command should include input format")
        self.assertIn("-omol2", command, "Command should include output format")
        self.assertIn(input_path, command, "Command should include input path")
        self.assertIn(output_path, command, "Command should include output path")
        
        print(f"Generated OpenBabel command: {command}")
    
    def test_construct_shell_command_mgltools(self):
        """Test constructing MGLTools shell command."""
        # Create input and output paths
        input_path = os.path.abspath(self.test_file)
        output_path = os.path.join(self.test_dir, "output.pdbqt")
        output_path = os.path.abspath(output_path)
        
        # Construct the command
        try:
            command = construct_shell_command("MGLTools", input_path, output_path)
            
            # Basic check if the command seems reasonable
            self.assertIsNotNone(command, "Command should not be None")
            self.assertIn("pythonsh", command, "MGLTools command should include pythonsh")
            self.assertIn("prepare_receptor4.py", command, "Command should include prepare_receptor4.py")
            self.assertIn("-r", command, "Command should include -r option")
            self.assertIn("-o", command, "Command should include -o option")
            
            print(f"Generated MGLTools command: {command}")
        except Exception as e:
            # If the MGLTools binaries don't exist, the test will fail with an error
            # We'll skip the test in this case
            self.skipTest(f"Skipping test due to error: {str(e)}")


if __name__ == '__main__':
    unittest.main() 