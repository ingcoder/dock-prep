import unittest
import os
import sys
import tempfile
import shutil

class TestFileOperations(unittest.TestCase):
    """Tests for file operations in PdbqtConverter."""
    
    def setUp(self):
        """Set up test environment before each test."""
        # Get the project root directory (2 levels up from this test file)
        project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        package_dir = "dock_prep"
        self.example_pdb = os.path.join(project_root, package_dir, "examples", "2pgh_original.pdb")
        self.test_dir = tempfile.mkdtemp()
        self.test_output_pdb = os.path.join(self.test_dir, "test_output.pdb")
    
    def tearDown(self):
        """Clean up after each test."""
        shutil.rmtree(self.test_dir)
    
    def test_extract_chains_manually(self):
        """Test manually extracting chains to a new PDB file."""
        # Read the example PDB file
        with open(self.example_pdb, 'r') as fin:
            lines = fin.readlines()
        
        # Find available chains
        chain_ids = set()
        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                chain_id = line[21]
                chain_ids.add(chain_id)
        
        if not chain_ids:
            self.skipTest("No chains found in example PDB")
        
        # Choose the first chain to extract
        target_chain = next(iter(sorted(chain_ids)))
        
        # Extract the chain
        with open(self.test_output_pdb, 'w') as fout:
            for line in lines:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    chain_id = line[21]
                    if chain_id == target_chain:
                        fout.write(line)
                elif line.startswith("TER") and len(line) > 21:
                    chain_id = line[21]
                    if chain_id == target_chain:
                        fout.write(line)
                else:
                    fout.write(line)
        
        # Verify the output file exists and has content
        self.assertTrue(os.path.exists(self.test_output_pdb), 
                       "Output PDB file should exist")
        self.assertGreater(os.path.getsize(self.test_output_pdb), 0, 
                         "Output PDB file should not be empty")
        
        # Verify the extracted file contains only the specified chain
        with open(self.test_output_pdb, 'r') as f:
            output_lines = f.readlines()
        
        extracted_chain_ids = set()
        for line in output_lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                chain_id = line[21]
                extracted_chain_ids.add(chain_id)
        
        self.assertEqual(extracted_chain_ids, {target_chain}, 
                        f"Extracted PDB should only contain chain {target_chain}")
        
        # Print statistics
        atoms_count = sum(1 for line in output_lines if line.startswith("ATOM"))
        hetatm_count = sum(1 for line in output_lines if line.startswith("HETATM"))
        print(f"Extracted chain {target_chain} with {atoms_count} ATOM and {hetatm_count} HETATM records")
    
    def test_modify_chain_ids(self):
        """Test modifying chain IDs in a PDB file."""
        # Create a copy of the example PDB for testing
        input_path = os.path.join(self.test_dir, "input.pdb")
        output_path = os.path.join(self.test_dir, "output.pdb")
        
        with open(self.example_pdb, 'r') as fin, open(input_path, 'w') as fout:
            lines = fin.readlines()[:200]  # Using a subset for testing
            fout.writelines(lines)
        
        # Identify original chain IDs
        chain_ids = set()
        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                chain_id = line[21]
                chain_ids.add(chain_id)
        
        if not chain_ids:
            self.skipTest("No chains found in example PDB")
        
        # Create a mapping from old to new chain IDs
        chain_mapping = {}
        for i, chain_id in enumerate(sorted(chain_ids)):
            # Use letters starting from X
            new_chain_id = chr(ord('X') + i)
            chain_mapping[chain_id] = new_chain_id
        
        # Modify chain IDs
        with open(input_path, 'r') as fin, open(output_path, 'w') as fout:
            for line in fin:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    chain_id = line[21]
                    if chain_id in chain_mapping:
                        new_chain_id = chain_mapping[chain_id]
                        line = line[:21] + new_chain_id + line[22:]
                fout.write(line)
        
        # Verify the modified chains
        with open(output_path, 'r') as f:
            modified_lines = f.readlines()
        
        modified_chain_ids = set()
        for line in modified_lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                chain_id = line[21]
                modified_chain_ids.add(chain_id)
        
        expected_modified_chains = set(chain_mapping.values())
        self.assertEqual(modified_chain_ids, expected_modified_chains, 
                       "Modified PDB should contain the mapped chain IDs")
        
        # Print the mapping
        print(f"Chain mapping: {chain_mapping}")
        for old_chain, new_chain in chain_mapping.items():
            old_count = sum(1 for line in lines if line.startswith(("ATOM", "HETATM")) and line[21] == old_chain)
            new_count = sum(1 for line in modified_lines if line.startswith(("ATOM", "HETATM")) and line[21] == new_chain)
            print(f"Chain {old_chain} -> {new_chain}: {old_count} atoms before, {new_count} atoms after")


if __name__ == '__main__':
    unittest.main() 