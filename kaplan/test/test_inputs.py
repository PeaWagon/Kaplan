"""Test the inputs module of Kaplan.

Note: each time the update_inputs(dict)
function is called, it resets the current
state of the borg to its defaults, which
ensures no previous values are being saved
between tests.

"""


import os

import vetee

from numpy.testing import assert_raises

from kaplan.inputs import Inputs, InputError


# directory for this test file
TEST_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'testfiles')


def test_inputs():
    test = Inputs()
    # create another instance of inputs
    test2 = Inputs()
    # change the method to something else
    test2.method = "dummy-method"
    # ensure other test instance receives the change
    assert test.method == "dummy-method"
    assert vars(test)
    # reset to default
    test._reset_to_defaults()

def test_inputs_update_inputs():
    """Test the update_inputs method of the Inputs class."""
    test = Inputs()
    # check that defaults have been re-established
    # from previous test
    assert test.method == "hf"
    # create minimum subset of inputs
    test_dict = {
        "struct_input": os.path.join(TEST_DIR, "1,3-butadiene.xyz"),
        "struct_type": "xyz",
        "charge": 0,
        "multip": 1,
        "method": "mp2"
    }
    # make sure no errors happen with inputs
    test.update_inputs(test_dict)
    # check method was updated
    assert test.method == "mp2"
    # create new instance of Inputs
    test2 = Inputs()
    # check that the method was updated for the new instance
    assert test2.struct_type == "xyz"
    # rest to defaults
    test._reset_to_defaults()
    # check reset was performed
    assert test2.method == "hf"
    # now try to input incorrect sets of inputs
    test_dict["struct_type"] = "not-avail"
    assert_raises(InputError, test.update_inputs, test_dict)
    # now try an invalid file name
    test_dict["struct_type"] = "xyz"
    test_dict["struct_input"] = "bloop"
    with assert_raises(FileNotFoundError):
        test.update_inputs(test_dict)
    # test with another type of input file
    test_dict["struct_input"] = os.path.join(TEST_DIR, "fopt-cis-1-chloropropene.com")
    # which won't work since the type is xyz and the input is com
    with assert_raises(InputError):
        test.update_inputs(test_dict)
    # now it should work
    test_dict["struct_type"] = "com"
    test.update_inputs(test_dict)
    
    # now try cases where incorrect program is given
    test_dict["prog"] = "not-avail"
    assert_raises(AssertionError, test.update_inputs, test_dict)
    test_dict["prog"] = "psi4"

    # cysteine
    test_dict["struct_type"] = "smiles"
    test_dict["struct_input"] = "C(C(C(=O)O)N)S"
    test.update_inputs(test_dict)
    assert test.num_atoms == 14
    assert "cysteine" in test.parser._synonyms

    # proline
    test_dict["struct_type"] = "name"
    test_dict["struct_input"] = "proline"
    test.update_inputs(test_dict)
    assert test.num_atoms == 17
    assert "proline" in test.parser._synonyms
    
    # threonine
    test_dict["struct_type"] = "cid"
    test_dict["struct_input"] = "6288"
    test.update_inputs(test_dict)
    assert test.num_atoms == 17
    assert "threonine" in test.parser._synonyms

    # test improper charge
    test_dict["charge"] = 0.3
    assert_raises(InputError, test.update_inputs, test_dict)
    # charge should work as float (but has to be whole number)
    test_dict["charge"] = 1.0
    test.update_inputs(test_dict)

    # test bad smiles string
    test_dict["charge"] = 1
    test_dict["struct_type"] = "smiles"
    test_dict["struct_input"] = "very-bad-smiles-string"
    assert_raises(vetee.structure.StructureError, test.update_inputs, test_dict)
    
    # test bad multiplicity
    test_dict["struct_input"] = "C=CCC=C"
    test_dict["multip"] = -3
    assert_raises(AssertionError, test.update_inputs, test_dict)

    # test molecule that is too small
    test_dict["struct_input"] = "hydrogen"
    test_dict["struct_type"] = "name"
    test_dict["multip"] = 1
    assert_raises(AssertionError, test.update_inputs, test_dict)
    test_dict["struct_input"] = "caffeine"

    # TODO: see why vetee is going into an infinite loop
    # when reading log files...
    #test_dict["struct_input"] = os.path.join(TEST_DIR, "1054_cc-pVDZ.log")
    #test.update_inputs(test_dict)

    # test bad ga inputs
    test_dict["num_slots"] = -1
    assert_raises(AssertionError, test.update_inputs, test_dict)
    test_dict["num_slots"] = 100

    test_dict["num_filled"] = 150
    assert_raises(AssertionError, test.update_inputs, test_dict)
    test_dict["num_filled"] = 20

    test_dict["num_mevs"] = -300
    assert_raises(AssertionError, test.update_inputs, test_dict)
    test_dict["num_mevs"] = 1000

    test_dict["num_swaps"] = 20
    assert_raises(AssertionError, test.update_inputs, test_dict)
    test_dict["num_swaps"] = 1

    test_dict["num_muts"] = 30
    assert_raises(AssertionError, test.update_inputs, test_dict)
    test_dict["num_muts"] = 3

    test_dict["num_geoms"] = -2
    assert_raises(AssertionError, test.update_inputs, test_dict)
    test_dict["num_geoms"] = 3

    test_dict["num_atoms"] = 2
    assert_raises(AssertionError, test.update_inputs, test_dict)
    test_dict["num_atoms"] = 24

    test_dict["fit_form"] = 1
    assert_raises(AssertionError, test.update_inputs, test_dict)
    test_dict["fit_form"] = 0

    test_dict["pmem_dist"] = 57
    assert_raises(AssertionError, test.update_inputs, test_dict)
    test_dict["pmem_dist"] = 5

    test_dict["coef_energy"] = -5
    assert_raises(AssertionError, test.update_inputs, test_dict)
    test_dict["coef_energy"] = 0.5

    test_dict["coef_rmsd"] = -5
    assert_raises(AssertionError, test.update_inputs, test_dict)
    test_dict["coef_rmsd"] = 0.5

    test_dict["t_size"] = 25
    assert_raises(AssertionError, test.update_inputs, test_dict)


test_inputs()
test_inputs_update_inputs()