#!/usr/bin/env python3
"""
Test Runner for PdbqtConverter

This script runs all the test files in the project.
"""

import unittest
import os
import sys

def main():
    # Discover and run all tests
    print("\n=== Running tests for PdbqtConverter ===")
    
    # Get the directory containing this script
    test_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Discover all test files in the current directory
    test_loader = unittest.TestLoader()
    test_suite = test_loader.discover(test_dir, pattern="test_*.py")
    
    # Run the tests
    test_runner = unittest.TextTestRunner(verbosity=2)
    result = test_runner.run(test_suite)
    
    # Print summary
    print("\n=== Test Summary ===")
    print(f"Tests run: {result.testsRun}")
    print(f"Errors: {len(result.errors)}")
    print(f"Failures: {len(result.failures)}")
    print(f"Skipped: {len(result.skipped)}")
    
    # Return non-zero exit code if any tests failed
    if result.errors or result.failures:
        print("\n❌ Some tests failed!")
        return 1
    else:
        print("\n✅ All tests passed!")
        return 0


if __name__ == "__main__":
    sys.exit(main()) 