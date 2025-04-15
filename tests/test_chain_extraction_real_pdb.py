import unittest
import os
import sys
import tempfile
import shutil

# Add the parent directory to the path so we can import the package
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from dock_prep.structure_handler import (
    save_clean_structure,
    load_structure_as_pdbfixer,
    map_pdbfixer_chains_to_original,
    get_chain_indices,
    get_chains_to_extract
)

"""
This module contains unit tests for PDB chain extraction functionality using real PDB files.
It tests the ability to identify, map, and extract protein chains (both ATOM and HETATM records)
based on various selection strategies including direct chain selection and distance-based
selection. Tests verify proper chain identification, mapping, and proximity-based extraction
with different cutoff distances.

Specifically, this test suite validates:
1. Chain mapping correctness between original PDB and PDBFixer representation
2. Retrieval of chain indices for single and multiple ATOM chains
3. Retrieval of chain indices for HETATM chains
4. Handling of non-existent chains and wrong record types
5. Direct inclusion of single and multiple chains
6. Using ATOM chains as reference for proximity-based extraction
7. Using HETATM chains as reference for proximity-based extraction
8. Effect of different distance cutoffs on chain extraction results
9. Chain extraction with various combinations of selection criteria

All tests utilize a real PDB file (2pgh_original.pdb) to ensure functionality works
with actual protein structures rather than simplified test cases.
"""

class TestChainExtractionRealPDB(unittest.TestCase):
    """Test chain extraction using real PDB files."""
    
    def setUp(self):
        """Set up test environment with real PDB file."""
        self.test_dir = tempfile.mkdtemp()
        project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        package_dir = "dock_prep"
        
        # Path to the example PDB file
        self.example_pdb = os.path.join(project_root, package_dir, "examples", "2pgh_original.pdb")
        
        # Create output files
        self.cleaned_pdb = os.path.join(self.test_dir, "cleaned.pdb")
        
        # First, clean the PDB file and get the original chain IDs
        self.original_chain_ids = save_clean_structure(self.example_pdb, self.cleaned_pdb, skip_hetatm=False, verbose=False)
        
        # Load the cleaned structure
        self.fixer = load_structure_as_pdbfixer(self.cleaned_pdb)
        
        # Create chain mapping
        self.chain_map = map_pdbfixer_chains_to_original(self.fixer, self.original_chain_ids)
        
        # Store actual chain IDs for verification
        self.atom_chains = []
        self.hetatm_chains = []
        for key in self.chain_map:
            if key[1] == "ATOM":
                self.atom_chains.append(key[0])
            elif key[1] == "HETATM":
                self.hetatm_chains.append(key[0])
        
        # print(f"PDB contains ATOM chains: {self.atom_chains}")
        #print(f"PDB contains HETATM chains: {self.hetatm_chains}")
    
    def tearDown(self):
        """Clean up temporary files."""
        shutil.rmtree(self.test_dir)
    
    def test_chain_map_correctness(self):
        """Test that the chain map correctly identifies chain types.
        
        VERIFIES:
        - Each chain map entry has the required fields (fixer_index, original_chain_id)
        - Chain IDs in keys match the original_chain_id values
        - All fixer_index values are valid integers
        - All record types are either ATOM or HETATM
        """
        print("\nTEST: Verifying chain map has correct structure and valid entries")
        # Verify that each chain in the map has the expected attributes
        for key, value in self.chain_map.items():
            chain_id, record_type = key
            self.assertIn("fixer_index", value, f"Chain map entry for {key} should have fixer_index")
            self.assertIn("original_chain_id", value, f"Chain map entry for {key} should have original_chain_id")
            
            # Verify that chain IDs match up
            self.assertEqual(chain_id, value["original_chain_id"], 
                             f"Chain ID in key {chain_id} should match original_chain_id {value['original_chain_id']}")
            
            # Verify that each fixer_index is a valid index
            self.assertIsInstance(value["fixer_index"], int, "fixer_index should be an integer")
            self.assertGreaterEqual(value["fixer_index"], 0, "fixer_index should be non-negative")
            
            # Verify record types are either ATOM or HETATM
            self.assertIn(record_type, ["ATOM", "HETATM"], "Record type must be either ATOM or HETATM")
        print(f"✅ RESULT: Successfully verified {len(self.chain_map)} chain map entries")
    
    def test_get_chain_indices_single_atom_chain(self):
        """Test retrieving chain index for a single ATOM chain.
        
        VERIFIES:
        - Can successfully retrieve indices for a single ATOM chain
        - Returned indices are valid integers
        """
        # Skip if no ATOM chains are available
        if not self.atom_chains:
            self.skipTest("No ATOM chains found in PDB file")
        
        # Get the first ATOM chain
        chain_id = self.atom_chains[0]
        
        print(f"\nTEST: Getting indices for single ATOM chain '{chain_id}'")
        # Get chain indices
        chain_indices = get_chain_indices(self.chain_map, [chain_id], "ATOM")
        self.assertTrue(len(chain_indices) > 0, f"Should find index for chain {chain_id}")
        
        # Verify that the chain indices are integers
        for idx in chain_indices:
            self.assertIsInstance(idx, int, "Chain index should be an integer")
        
        print(f"✅ RESULT: Found {len(chain_indices)} indices for ATOM chain '{chain_id}': {chain_indices}")
    
    def test_get_chain_indices_multiple_atom_chains(self):
        """Test retrieving chain indices for multiple ATOM chains.
        
        VERIFIES:
        - Can successfully retrieve indices for multiple ATOM chains
        - Correct number of indices are returned
        """
        # Skip if we don't have at least 2 ATOM chains
        if len(self.atom_chains) < 2:
            self.skipTest("Need at least 2 ATOM chains for this test")
        
        # Get the first two ATOM chains
        chain_ids = self.atom_chains[:2]
        
        print(f"\nTEST: Getting indices for multiple ATOM chains {chain_ids}")
        # Get chain indices
        chain_indices = get_chain_indices(self.chain_map, chain_ids, "ATOM")
        self.assertEqual(len(chain_indices), 2, f"Should find exactly 2 indices for chains {chain_ids}")
        
        print(f"✅ RESULT: Found indices for chains {chain_ids}: {chain_indices}")
    
    def test_get_chain_indices_hetatm_chain(self):
        """Test retrieving chain index for a HETATM chain.
        
        VERIFIES:
        - Can successfully retrieve indices for HETATM chains
        """
        # Skip if no HETATM chains are available
        if not self.hetatm_chains:
            self.skipTest("No HETATM chains found in PDB file")
        
        # Get the first HETATM chain
        chain_id = self.hetatm_chains[0]
        
        print(f"\nTEST: Getting indices for HETATM chain '{chain_id}'")
        # Get chain indices
        chain_indices = get_chain_indices(self.chain_map, [chain_id], "HETATM")
        self.assertTrue(len(chain_indices) > 0, f"Should find index for chain {chain_id}")
        
        print(f"✅ RESULT: Found {len(chain_indices)} indices for HETATM chain '{chain_id}': {chain_indices}")
    
    def test_get_chain_indices_nonexistent_chain(self):
        """Test retrieving chain index for a non-existent chain.
        
        VERIFIES:
        - Returns empty list when requesting a chain that doesn't exist
        """
        # Try a non-existent chain ID
        chain_id = "Z"  # Assuming Z doesn't exist in the PDB
        
        print(f"\nTEST: Getting indices for non-existent chain '{chain_id}'")
        # Get chain indices
        chain_indices = get_chain_indices(self.chain_map, [chain_id], "ATOM")
        self.assertEqual(chain_indices, [], f"Should return empty list for non-existent chain {chain_id}")
        
        print(f"✅ RESULT: Correctly returned empty list for non-existent chain '{chain_id}'")
    
    def test_get_chain_indices_wrong_record_type(self):
        """Test retrieving chain index with incorrect record type.
        
        VERIFIES:
        - Returns empty list when requesting an ATOM chain as HETATM (or vice versa)
        - Correctly handles chains that exist as both types
        """
        # Skip if no ATOM chains are available
        if not self.atom_chains:
            self.skipTest("No ATOM chains found in PDB file")
        
        # Get the first ATOM chain but request it as HETATM
        chain_id = self.atom_chains[0]
        
        print(f"\nTEST: Getting indices for ATOM chain '{chain_id}' but requested as HETATM")
        # Get chain indices with the wrong record type
        chain_indices = get_chain_indices(self.chain_map, [chain_id], "HETATM")
        
        # Check if this chain also exists as a HETATM chain
        is_also_hetatm = chain_id in self.hetatm_chains
        
        if is_also_hetatm:
            self.assertTrue(len(chain_indices) > 0, 
                          f"Chain {chain_id} exists as both ATOM and HETATM, should find an index")
            print(f"✅ RESULT: Chain '{chain_id}' exists as both ATOM and HETATM, found indices: {chain_indices}")
        else:
            self.assertEqual(chain_indices, [], 
                           f"Should return empty list for chain {chain_id} with wrong record type")
            print(f"✅ RESULT: Correctly returned empty list for ATOM chain '{chain_id}' when requested as HETATM")
    
    def test_include_single_chain(self):
        """Test including a single chain directly.
        
        VERIFIES:
        - Can extract a single chain by direct inclusion
        - Results have the correct format (chain_id, record_type)
        - Record types are valid (ATOM or HETATM)
        """
        # Skip if no ATOM chains are available
        if not self.atom_chains:
            self.skipTest("No ATOM chains found in PDB file")
        
        # Get the first ATOM chain
        chain_id = self.atom_chains[0]
        
        print(f"\nTEST: Including single ATOM chain '{chain_id}' directly")
        # Get chain indices
        chain_indices = get_chain_indices(self.chain_map, [chain_id], "ATOM")
        self.assertTrue(len(chain_indices) > 0, f"Should find index for chain {chain_id}")
        
        # Get chains to extract
        chains_to_extract = get_chains_to_extract(self.chain_map, chain_indices, self.fixer, 5.0)
        self.assertTrue(len(chains_to_extract) > 0, "Should extract at least one chain")
        
        # Verify the result format
        for chain_tuple in chains_to_extract:
            self.assertEqual(len(chain_tuple), 2, "Each chain should be a tuple of (chain_id, record_type)")
            self.assertIn(chain_tuple[1], ["ATOM", "HETATM"], "Record type should be ATOM or HETATM")
        
        # Log chains for reference
        chain_ids_in_result = [f"{chain[0]} ({chain[1]})" for chain in chains_to_extract]
        print(f"✅ RESULT: Using INCLUDE_CHAIN_IDS='{chain_id}', extracted chains: {chain_ids_in_result}")
    
    def test_include_multiple_chains(self):
        """Test including multiple chains directly.
        
        VERIFIES:
        - Can extract multiple chains by direct inclusion
        - At least one target chain appears in the results
        """
        # Skip if we don't have at least 2 ATOM chains
        if len(self.atom_chains) < 2:
            self.skipTest("Need at least 2 ATOM chains for this test")
        
        # Get the first two ATOM chains
        chain_ids = self.atom_chains[:2]
        
        print(f"\nTEST: Including multiple ATOM chains {chain_ids} directly")
        # Get chain indices
        chain_indices = get_chain_indices(self.chain_map, chain_ids, "ATOM")
        self.assertEqual(len(chain_indices), 2, f"Should find indices for chains {chain_ids}")
        
        # Get chains to extract
        chains_to_extract = get_chains_to_extract(self.chain_map, chain_indices, self.fixer, 5.0)
        self.assertTrue(len(chains_to_extract) > 0, "Should extract at least one chain")
        
        # Verify the result contains at least one of our target chains
        chain_ids_in_result = [chain[0] for chain in chains_to_extract]
        self.assertTrue(any(chain_id in chain_ids_in_result for chain_id in chain_ids),
                      f"At least one of {chain_ids} should be in the result: {chain_ids_in_result}")
        
        # Log chains for reference
        chain_ids_in_result = [f"{chain[0]} ({chain[1]})" for chain in chains_to_extract]
        print(f"✅ RESULT: Using INCLUDE_CHAIN_IDS={','.join(chain_ids)}, extracted chains: {chain_ids_in_result}")
    
    def test_reference_atom_chain(self):
        """Test using ATOM chain as reference.
        
        VERIFIES:
        - Can use a single ATOM chain as reference for extraction
        - At least one chain is extracted (which may be the reference chain)
        """
        # Skip if no ATOM chains are available
        if not self.atom_chains:
            self.skipTest("No ATOM chains found in PDB file")
        
        # Get the first ATOM chain
        chain_id = self.atom_chains[0]
        
        print(f"\nTEST: Using ATOM chain '{chain_id}' as reference with 5.0Å cutoff")
        # Get chain indices
        chain_indices = get_chain_indices(self.chain_map, [chain_id], "ATOM")
        self.assertTrue(len(chain_indices) > 0, f"Should find index for chain {chain_id}")
        
        # Get chains to extract using ATOM as reference
        chains_to_extract = get_chains_to_extract(self.chain_map, chain_indices, self.fixer, 5.0, filter_type="ATOM")
        
        # Verify the result (should have at least the reference chain or other chains within the cutoff)
        self.assertTrue(len(chains_to_extract) > 0, "Should extract at least one chain")
        
        # Log results for reference
        chain_ids_in_result = [f"{chain[0]} ({chain[1]})" for chain in chains_to_extract]
        print(f"✅ RESULT: Using REFERENCE_ATOM_CHAINS='{chain_id}', extracted chains: {chain_ids_in_result}")
    
    def test_reference_atom_chains_multiple(self):
        """Test using multiple ATOM chains as reference.
        
        VERIFIES:
        - Can use multiple ATOM chains as reference for extraction
        - At least one chain is extracted
        """
        # Skip if we don't have at least 2 ATOM chains
        if len(self.atom_chains) < 2:
            self.skipTest("Need at least 2 ATOM chains for this test")
        
        # Get the first two ATOM chains
        chain_ids = self.atom_chains[:2]
        
        # Get chain indices
        chain_indices = get_chain_indices(self.chain_map, chain_ids, "ATOM")
        self.assertEqual(len(chain_indices), 2, f"Should find indices for chains {chain_ids}")
        
        # Get chains to extract using ATOM chains as reference
        chains_to_extract = get_chains_to_extract(self.chain_map, chain_indices, self.fixer, 5.0, filter_type="ATOM")
        
        # Verify the result (should have chains within the cutoff)
        self.assertTrue(len(chains_to_extract) > 0, "Should extract at least one chain")
        
        # Log results for reference
        chain_ids_in_result = [f"{chain[0]} ({chain[1]})" for chain in chains_to_extract]
        print(f"✅ RESULT: Using REFERENCE_ATOM_CHAINS={','.join(chain_ids)}, found chains: {chain_ids_in_result}")
    
    def test_reference_hetatm_chain(self):
        """Test using HETATM chain as reference.
        
        VERIFIES:
        - Can use a HETATM chain as reference for extraction
        - Either finds nearby ATOM chains or returns empty if none are close enough
        """
        # Skip if no HETATM chains are available
        if not self.hetatm_chains:
            self.skipTest("No HETATM chains found in PDB file")
        
        # Get the first HETATM chain
        chain_id = self.hetatm_chains[0]
        
        # Get chain indices
        chain_indices = get_chain_indices(self.chain_map, [chain_id], "HETATM")
        self.assertTrue(len(chain_indices) > 0, f"Should find index for chain {chain_id}")
        
        # Get chains to extract using HETATM as reference
        chains_to_extract = get_chains_to_extract(self.chain_map, chain_indices, self.fixer, 5.0, filter_type="HETATM")
        
        # Verify the result (should have at least one ATOM chain within the cutoff)
        chain_types_in_result = [chain[1] for chain in chains_to_extract]
        has_atom_chains = "ATOM" in chain_types_in_result
        
        # In a real PDB, we expect to find at least one ATOM chain near HETATM chains
        self.assertTrue(has_atom_chains or len(chains_to_extract) == 0, 
                      "Should have either ATOM chains or no chains at all")
        
        # Log results for reference
        chain_ids_in_result = [f"{chain[0]} ({chain[1]})" for chain in chains_to_extract]
        print(f"✅ RESULT: Using REFERENCE_HETATM_CHAINS={chain_id}, found chains: {chain_ids_in_result}")
    
    def test_reference_hetatm_chains_multiple(self):
        """Test using multiple HETATM chains as reference.
        
        VERIFIES:
        - Can use multiple HETATM chains as reference for extraction
        - Either finds nearby ATOM chains or returns empty if none are close enough
        """
        # Skip if we don't have at least 2 HETATM chains
        if len(self.hetatm_chains) < 2:
            self.skipTest("Need at least 2 HETATM chains for this test")
        
        # Get the first two HETATM chains
        chain_ids = self.hetatm_chains[:2]
        
        # Get chain indices
        chain_indices = get_chain_indices(self.chain_map, chain_ids, "HETATM")
        self.assertEqual(len(chain_indices), 2, f"Should find indices for chains {chain_ids}")
        
        # Get chains to extract using HETATM chains as reference
        chains_to_extract = get_chains_to_extract(self.chain_map, chain_indices, self.fixer, 5.0, filter_type="HETATM")
        
        # Verify the result (should have at least one ATOM chain within the cutoff)
        chain_types_in_result = [chain[1] for chain in chains_to_extract]
        has_atom_chains = "ATOM" in chain_types_in_result
        
        # In a real PDB, we expect to find at least one ATOM chain near HETATM chains
        self.assertTrue(has_atom_chains or len(chains_to_extract) == 0, 
                      "Should have either ATOM chains or no chains at all")
        
        # Log results for reference
        chain_ids_in_result = [f"{chain[0]} ({chain[1]})" for chain in chains_to_extract]
        print(f"✅ RESULT: Using REFERENCE_HETATM_CHAINS={','.join(chain_ids)}, found chains: {chain_ids_in_result}")
    
    def test_cutoff_distance_effect(self):
        """Test the effect of different cutoff distances.
        
        VERIFIES:
        - Larger cutoff distance finds at least as many chains as smaller cutoff
        - Distance cutoff directly affects number of chains extracted
        """
        # Skip if no HETATM chains are available
        if not self.hetatm_chains:
            self.skipTest("No HETATM chains found in PDB file")
        
        # Get the first HETATM chain
        chain_id = self.hetatm_chains[0]
        
        print(f"\nTEST: Testing different cutoff distances using HETATM chain '{chain_id}' as reference")
        # Get chain indices
        chain_indices = get_chain_indices(self.chain_map, [chain_id], "HETATM")
        
        # Test with 3.0 Å cutoff (smaller)
        small_cutoff = 3.0
        print(f"Testing with small cutoff distance: {small_cutoff}Å")
        chains_small_cutoff = get_chains_to_extract(self.chain_map, chain_indices, self.fixer, small_cutoff, filter_type="HETATM")
        
        # Test with 10.0 Å cutoff (larger)
        large_cutoff = 10.0
        print(f"Testing with large cutoff distance: {large_cutoff}Å")
        chains_large_cutoff = get_chains_to_extract(self.chain_map, chain_indices, self.fixer, large_cutoff, filter_type="HETATM")
        
        # The larger cutoff should find at least as many chains as the smaller cutoff
        self.assertTrue(len(chains_large_cutoff) >= len(chains_small_cutoff),
                       f"Larger cutoff {large_cutoff}Å should find at least as many chains as smaller cutoff {small_cutoff}Å")
        
        # Log results
        small_cutoff_chains = [f"{chain[0]} ({chain[1]})" for chain in chains_small_cutoff]
        large_cutoff_chains = [f"{chain[0]} ({chain[1]})" for chain in chains_large_cutoff]
        print(f"✅ RESULT: With {small_cutoff}Å cutoff: Found {len(small_cutoff_chains)} chains: {small_cutoff_chains}")
        print(f"✅ RESULT: With {large_cutoff}Å cutoff: Found {len(large_cutoff_chains)} chains: {large_cutoff_chains}")
        print(f"CONCLUSION: Larger cutoff found {len(large_cutoff_chains) - len(small_cutoff_chains)} more chains")

if __name__ == "__main__":
    unittest.main() 