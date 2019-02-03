"""This test runs the main Kaplan program
with some dummy input files."""

import os

from kaplan.gac import run_kaplan

# directory for this test file
TEST_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'testfiles')


def test_run_kaplan():
    """Test run_kaplan function from gac module."""
    test_ga = os.path.join(TEST_DIR, "example_ga_input_file.txt")
    test_mol = os.path.join(TEST_DIR, "example_mol_input_file.txt")
    run_kaplan(test_ga, test_mol)
    test_ga2 = os.path.join(TEST_DIR, "example2_ga_input_file.txt")
    test_mol2 = os.path.join(TEST_DIR, "example2_mol_input_file.txt")
    run_kaplan(test_ga2, test_mol2)
