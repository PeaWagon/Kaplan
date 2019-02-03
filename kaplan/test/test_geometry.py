"""Test geometry module from Kaplan."""

import os

# from numpy.testing import assert_raises

from kaplan.geometry import generate_parser
# get_zmatrix_template, update_zmatrix, zmatrix_to_xyz


# directory for this test file
TEST_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'testfiles')


def test_generate_parser():
    """Test generate_parser function from geometry module."""
    # ensure valid inputs can be used for parser object
    mol_input_dict = {"struct_type": "name", "struct_input": "water",
                      "qcm": "hf", "basis": "sto-3g", "charge": 0, "multip": 1}
    generate_parser(mol_input_dict)
    mol_input_dict["struct_type"] = "cid"
    mol_input_dict["struct_input"] = "8005"
    generate_parser(mol_input_dict)
    mol_input_dict["struct_type"] = "smiles"
    mol_input_dict["struct_input"] = "c1=cc=cc=c1"
    generate_parser(mol_input_dict)
    mol_input_dict["struct_type"] = "xyz"
    mol_input_dict["struct_input"] = os.path.join(TEST_DIR, "1,3-butadiene.xyz")
    generate_parser(mol_input_dict)
    mol_input_dict["struct_type"] = "com"
    mol_input_dict["struct_input"] = os.path.join(TEST_DIR, "fopt-cis-1-chloropropene.com")
    generate_parser(mol_input_dict)
    mol_input_dict["struct_type"] = "glog"
    mol_input_dict["struct_input"] = os.path.join(TEST_DIR, "fopt-cis-1-chloro-2-fluoroethene.log")
    generate_parser(mol_input_dict)


def test_get_zmatrix_template():
    """Test get_zmatrix_template function from geometry module."""


def test_update_zmatrix():
    """Test update_zmatrix function from geometry module."""


def test_zmatrix_to_xyz():
    """Test zmatrix_to_xyz function from geometry module."""
