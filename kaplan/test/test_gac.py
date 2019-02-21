"""This test runs the main Kaplan program
with some dummy input files."""

import os

from kaplan.gac import run_kaplan

# directory for this test file
TEST_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'testfiles')


def test_run_kaplan():
    """Test run_kaplan function from gac module."""
    test1 = {
        "struct_input": os.path.join(TEST_DIR, "caffeine.xyz"),
        "struct_type": "xyz",
        "charge": 0,
        "multip": 1
    }
    run_kaplan(test1)
