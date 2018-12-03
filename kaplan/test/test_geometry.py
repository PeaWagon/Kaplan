
import os

from numpy.testing import assert_raises

from kaplan.geometry import generate_parser, get_zmatrix_template, update_zmatrix, zmatrix_to_xyz

# directory for this test file
test_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'testfiles')

def test_generate_parser():
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
    mol_input_dict["struct_input"] = os.path.join(test_dir, "1,3-butadiene.xyz")
    generate_parser(mol_input_dict)
    mol_input_dict["struct_type"] = "com"
    mol_input_dict["struct_input"] = os.path.join(test_dir, "fopt-cis-1-chloropropene.com")
    generate_parser(mol_input_dict)
    mol_input_dict["struct_type"] = "glog"
    mol_input_dict["struct_input"] = os.path.join(test_dir, "fopt-cis-1-chloro-2-fluoroethene.log")
    generate_parser(mol_input_dict)
    try:
        pass
    except Exception:
        pass

def test_get_zmatrix_template():
    pass

def test_update_zmatrix():
    pass

def test_zmatrix_to_xyz():
    pass

