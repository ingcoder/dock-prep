import unittest
import os
import sys
import tempfile
import shutil
from dock_prep.structure_handler import (
    load_clean_structure,
    save_structure_to_pdb,
    extract_chains_to_pdb,
    restore_original_chain_ids,
    get_missing_residues_by_chain
)

class TestStructureHandler(unittest.TestCase):
    """Test cases for dock-prep structure_handler module."""
    
    def setUp(self):
        """Set up test environment before each test."""
        project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        package_dir = "dock_prep"
        self.example_pdb = os.path.join(project_root, package_dir, "examples", "2pgh_original.pdb")
        self.test_dir = tempfile.mkdtemp()
        self.test_output_pdb = os.path.join(self.test_dir, "test_output.pdb")
        
    def tearDown(self):
        """Clean up after each test."""
        shutil.rmtree(self.test_dir)
    
    def test_load_clean_structure(self):
        """Test loading and cleaning a PDB structure."""
        # Test loading with no_hetatm=True
        fixer = load_clean_structure(self.example_pdb, no_hetatm=True, verbose=False)
        self.assertIsNotNone(fixer, "PDBFixer object should not be None")
        
        # Verify that the fixer object has a topology attribute
        self.assertTrue(hasattr(fixer, 'topology'), "PDBFixer object has notopology attribute")
        
        # Verify chains in topology
        chain_count = sum(1 for _ in fixer.topology.chains())
        self.assertGreater(chain_count, 0, "PDBFixer topology has one chain")
    
    def test_extract_chains_to_pdb(self):
        """Test extracting specific chains from a PDB file."""
        # First, get the chain IDs from the example PDB
        with open(self.example_pdb, 'r') as f:
            lines = f.readlines()
        
        chain_ids = set()
        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                chain_id = line[21]
                chain_ids.add(chain_id)
        
        if not chain_ids:
            self.skipTest("No chains found in example PDB")
        
        # Extract the first chain
        first_chain_id = next(iter(chain_ids))
        target_chains = [first_chain_id]
        
        extract_chains_to_pdb(self.example_pdb, self.test_output_pdb, target_chains)
        
        # Verify the extracted file contains only the specified chain
        with open(self.test_output_pdb, 'r') as f:
            lines = f.readlines()
        
        extracted_chain_ids = set()
        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                chain_id = line[21]
                extracted_chain_ids.add(chain_id)
        
        self.assertEqual(extracted_chain_ids, set(target_chains), 
                        f"Extracted PDB should only contain chains {target_chains}")
    
    def test_save_structure_to_pdb(self):
        """Test saving a PDBFixer structure to a PDB file."""
        fixer = load_clean_structure(self.example_pdb, verbose=False)
        save_structure_to_pdb(fixer, self.test_output_pdb, verbose=False)
        
        # Verify the file exists and has content
        self.assertTrue(os.path.exists(self.test_output_pdb), 
                      "Output PDB file should exist")
        self.assertGreater(os.path.getsize(self.test_output_pdb), 0, 
                         "Output PDB file should not be empty")
    
    def test_restore_original_chain_ids(self):
        """Test restoring original chain IDs in a processed PDB file."""
        # First, create a sample PDB with modified chain IDs
        input_path = os.path.join(self.test_dir, "modified_chains.pdb")
        output_path = os.path.join(self.test_dir, "restored_chains.pdb")
        
        # Get the first few lines from the example PDB and modify chain IDs
        with open(self.example_pdb, 'r') as f:
            lines = f.readlines()[:100]  # Get a subset for testing
        
        original_chain_ids = set()
        for i, line in enumerate(lines):
            if line.startswith("ATOM") or line.startswith("HETATM"):
                original_chain_id = line[21]
                original_chain_ids.add(original_chain_id)
                # Change chain ID to 'X' for testing
                lines[i] = line[:21] + 'X' + line[22:]
        
        with open(input_path, 'w') as f:
            f.writelines(lines)
        
        # Restore original chain IDs
        if original_chain_ids:
            target_chains = list(original_chain_ids)
            restore_original_chain_ids(input_path, output_path, target_chains, verbose=False)
            
            # Verify the restored file has the original chain IDs
            with open(output_path, 'r') as f:
                restored_lines = f.readlines()
            
            restored_chain_ids = set()
            for line in restored_lines:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    chain_id = line[21]
                    restored_chain_ids.add(chain_id)
            
            self.assertEqual(restored_chain_ids, original_chain_ids, 
                           "Restored PDB should have the original chain IDs")


if __name__ == '__main__':
    unittest.main() 