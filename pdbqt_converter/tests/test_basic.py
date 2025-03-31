import unittest
import os
import sys

class TestBasic(unittest.TestCase):
    """Basic tests for PdbqtConverter."""
    
    def setUp(self):
        """Set up test environment before each test."""
        self.example_pdb = os.path.join("examples", "2pgh_original.pdb")
    
    def test_example_file_exists(self):
        """Test that the example PDB file exists."""
        self.assertTrue(os.path.exists(self.example_pdb), 
                       f"Example PDB file {self.example_pdb} should exist")
    
    def test_file_content(self):
        """Test that the example PDB file has expected content."""
        with open(self.example_pdb, 'r') as f:
            lines = f.readlines()
        
        # Check if file has content
        self.assertGreater(len(lines), 0, "Example PDB file should not be empty")
        
        # Check if file has ATOM records
        atom_lines = [line for line in lines if line.startswith("ATOM")]
        self.assertGreater(len(atom_lines), 0, "Example PDB file should have ATOM records")
    
    def test_file_parsing(self):
        """Test basic file parsing of PDB format."""
        with open(self.example_pdb, 'r') as f:
            lines = f.readlines()
        
        # Extract chain IDs from ATOM records
        chain_ids = set()
        for line in lines:
            if line.startswith("ATOM"):
                chain_id = line[21]
                chain_ids.add(chain_id)
        
        # Check if chains were found
        self.assertGreater(len(chain_ids), 0, "Example PDB file should have at least one chain")
        
        # Print the found chains (helpful for debugging)
        print(f"Found chains: {', '.join(sorted(chain_ids))}")
    
    def test_simple_extract_chains(self):
        """Test extracting chains without using the structure_handler module."""
        with open(self.example_pdb, 'r') as f:
            lines = f.readlines()
        
        # Get all chains
        chain_ids = set()
        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                chain_id = line[21]
                chain_ids.add(chain_id)
        
        if not chain_ids:
            self.skipTest("No chains found in example PDB")
        
        # Extract the first chain
        first_chain_id = next(iter(sorted(chain_ids)))
        
        # Filter lines for the first chain
        filtered_lines = [line for line in lines if not (line.startswith("ATOM") or line.startswith("HETATM")) 
                          or line[21] == first_chain_id]
        
        # Count ATOM records for the specific chain
        atom_count = sum(1 for line in filtered_lines if line.startswith("ATOM") and line[21] == first_chain_id)
        
        # Check if extraction worked
        self.assertGreater(atom_count, 0, f"Should find ATOM records for chain {first_chain_id}")
        print(f"Found {atom_count} ATOM records for chain {first_chain_id}")


if __name__ == '__main__':
    unittest.main() 