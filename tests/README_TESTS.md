# Test Suite for PdbqtConverter

This directory contains several test files to verify the functionality of the dock-prep project. The tests range from basic file operations to more complex tests requiring external dependencies.

## Available Tests

1. **test_basic.py** - Basic tests that verify file existence and simple parsing
2. **test_file_operations.py** - Tests for file operations like chain extraction and ID modification
3. **test_subprocess_handler.py** - Tests for the subprocess_handler module
4. **test_structure_handler.py** - Tests for the structure_handler module (requires dependencies)
5. **run_tests.py** - A script to run all the tests at once

## Running Tests

### Basic Tests (No Dependencies Required)

To run the basic tests that don't require dependencies:

```bash
python3 -m unittest test_basic.py
python3 -m unittest test_file_operations.py
python3 -m unittest test_subprocess_handler.py
```

### Full Test Suite (Dependencies Required)

To run the full test suite, first install the required dependencies:

```bash
# Using pip for Python packages
pip install numpy biopython

# Using conda for pdbfixer and openmm
conda install -c conda-forge pdbfixer openmm
```

Then run all tests with:

```bash
python3 run_tests.py
```

Or you can run individual test files with:

```bash
python3 -m unittest test_structure_handler.py
```

## Adding New Tests

To add new tests:

1. Create a new file named `test_*.py` (the prefix `test_` is important)
2. Write your tests using the Python's unittest framework
3. The test runner will automatically discover and run your tests

## Best Practices

- Always use `setUp()` and `tearDown()` methods to set up and clean up test environments
- Use temporary directories for test files to avoid leaving residual files
- Use descriptive test method names (they should start with `test_`)
- Add docstrings to test methods to explain what they are testing

## Troubleshooting

If you encounter errors about missing dependencies, make sure you have installed all the required packages as mentioned above.

If you see errors related to file paths, check that the test is looking for files in the correct location. The tests assume they are run from the project root directory.

## Continuous Integration

In a future update, we plan to integrate these tests with a CI/CD pipeline to automatically run tests on every commit. 