"""
These tests are for the mol_input module. The mol_input
module has two functions: read_mol_input and verify_mol_input.
If these tests do not pass, then it is possible that
the program's constraints have been altered or the program
uses another program to construct the Parser object (or
vetee has been updated).

This program uses some testfiles located in the
testfiles folder of kaplan. Specifically, in
"../Kaplan/kaplan/test/testfiles", we have:
"example_mol_input_file.txt" (a good input file)
"example2_mol_input_file.txt" (a good input file with
newlines between the inputs)
"bad1_mol_input_file.txt" (too many arguments)
"bad2_mol_input_file.txt" (too few arguments)
"bad3_mol_input_file.txt" (incorrect spelling for argument)

"""

import os

from numpy.testing import assert_raises

from kaplan.mol_input import read_mol_input, verify_mol_input

# directory for this test file
test_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'testfiles')

def test_read_mol_input():
    # good input
    read_mol_input(os.path.join(test_dir, "example_mol_input_file.txt"))
    # good input with extra spaces
    read_mol_input(os.path.join(test_dir, "example2_mol_input_file.txt"))
    # no such file error
    assert_raises(FileNotFoundError, read_mol_input, 'no-such-file')
    # qcm appears twice
    assert_raises(ValueError, read_mol_input, os.path.join(test_dir, "bad1_mol_input_file.txt"))
    # missing struct type
    assert_raises(ValueError, read_mol_input, os.path.join(test_dir, "bad2_mol_input_file.txt"))

def test_verify_mol_input():
    # good input with extra spaces
    mol_input_dict = read_mol_input(os.path.join(test_dir, "example2_mol_input_file.txt"))
    print(mol_input_dict)
    verify_mol_input(mol_input_dict)

    # check the good file
    # {'qcm': 'hf', 'basis': 'sto-3g', 'struct_input': 'c=cc=c',
    #  'struct_type': 'smiles', 'prog': 'psi4', 'charge': '0', 'multip': '1'}
    mol_input_dict = read_mol_input(os.path.join(test_dir, "example_mol_input_file.txt"))
    verify_mol_input(mol_input_dict)

    # struct input is spelt wrong
    mol_input_dict2 = read_mol_input(os.path.join(test_dir, "bad3_mol_input_file.txt"))
    assert_raises(ValueError, verify_mol_input, mol_input_dict2)

    mol_input_dict["qcm"] = "not-a-method"
    assert_raises(ValueError, verify_mol_input, mol_input_dict)
    mol_input_dict["qcm"] = "hf"

    mol_input_dict["basis"] = "not-a-basis"
    assert_raises(ValueError, verify_mol_input, mol_input_dict)
    mol_input_dict["basis"] = "sto-3g"

    mol_input_dict["struct_input"] = "very-bad-smiles-string"
    assert_raises(ValueError, verify_mol_input, mol_input_dict)
    mol_input_dict["struct_input"] = "c=cc=c"

    mol_input_dict["struct_type"] = "not-an-option"
    assert_raises(AssertionError, verify_mol_input, mol_input_dict)
    mol_input_dict["struct_type"] = "smiles"

    mol_input_dict["prog"] = "unavailable-prog"
    assert_raises(AssertionError, verify_mol_input, mol_input_dict)
    mol_input_dict["prog"] = "psi4"

    mol_input_dict["charge"] = "0.34"
    assert_raises(ValueError, verify_mol_input, mol_input_dict)
    mol_input_dict["charge"] = "0"

    mol_input_dict["multip"] = "-2"
    assert_raises(AssertionError, verify_mol_input, mol_input_dict)

if __name__ == "__main__":
    test_read_mol_input()
    test_verify_mol_input()
