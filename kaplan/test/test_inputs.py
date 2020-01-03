"""Test the inputs module of Kaplan.

Note: each time the update_inputs(dict)
function is called, it resets the current
state of the borg to its defaults, which
ensures no previous values are being saved
between tests.

"""


import os
import pickle
import shutil

from copy import deepcopy

from numpy.testing import assert_raises

from kaplan.inputs import Inputs, InputError, read_input, get_latest_job
from kaplan.tools import TEST_DIR
from kaplan.control import run_kaplan
from kaplan.geometry import GeometryError


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
    out_dir = os.path.join(TEST_DIR, "test_update_inputs")

    # remove the test directory if it already exists
    if os.path.isdir(out_dir):
        shutil.rmtree(out_dir, ignore_errors=True)
    os.mkdir(out_dir)

    test = Inputs()
    test._reset_to_defaults()
    # check that defaults have been re-established
    # from previous test
    assert test.method is None
    # create minimum subset of inputs
    test_dict = {
        "struct_input": os.path.join(TEST_DIR, "1,3-butadiene.xyz"),
        "struct_type": "xyz",
        "charge": 0,
        "multip": 1,
        "method": "mp2",
        "output_dir": out_dir,
    }

    # should get an error since mp2 is not a forcefield
    # the default is to use openbabel
    with assert_raises(InputError):
        test.update_inputs(test_dict)
    test_dict["prog"] = "psi4"
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
    assert test2.method is None
    # now try to input incorrect sets of inputs
    test_dict["struct_type"] = "not-avail"
    assert_raises(InputError, test.update_inputs, test_dict)
    # now try an invalid file name
    test_dict["struct_type"] = "xyz"
    test_dict["struct_input"] = "bloop"
    with assert_raises(GeometryError):
        test.update_inputs(test_dict)
    # test with another type of input file
    test_dict["struct_input"] = os.path.join(TEST_DIR, "fopt-cis-1-chloropropene.com")
    # which won't work since the type is xyz and the input is com
    with assert_raises(GeometryError):
        test.update_inputs(test_dict)
    test_dict["struct_type"] = "com"
    with assert_raises(GeometryError):
        test.update_inputs(test_dict)

    # now try cases where incorrect program is given
    test_dict["prog"] = "not-avail"
    assert_raises(GeometryError, test.update_inputs, test_dict)
    test_dict["prog"] = "psi4"

    # cysteine
    test_dict["struct_type"] = "smiles"
    test_dict["struct_input"] = "C(C(C(=O)O)N)S"
    test.update_inputs(test_dict)
    assert len(test.coords) == 14

    # proline
    test_dict["struct_type"] = "name"
    test_dict["struct_input"] = "proline"
    test.update_inputs(test_dict)
    assert len(test.coords) == 17

    # threonine
    test_dict["struct_type"] = "cid"
    test_dict["struct_input"] = "6288"
    test.update_inputs(test_dict)
    assert len(test.coords) == 17

    # test improper charge
    test_dict["charge"] = 0.3
    assert_raises(InputError, test.update_inputs, test_dict)
    test_dict["charge"] = 1.0
    assert_raises(InputError, test.update_inputs, test_dict)

    # test bad smiles string
    test_dict["charge"] = 1
    test_dict["struct_type"] = "smiles"
    test_dict["struct_input"] = "very-bad-smiles-string"
    assert_raises(GeometryError, test.update_inputs, test_dict)

    # test bad multiplicity
    test_dict["struct_input"] = "C=CCC=C"
    test_dict["multip"] = -3
    assert_raises(InputError, test.update_inputs, test_dict)

    # test molecule that is too small
    test_dict["struct_input"] = "hydrogen"
    test_dict["struct_type"] = "name"
    test_dict["multip"] = 1
    assert_raises(InputError, test.update_inputs, test_dict)
    test_dict["struct_input"] = "caffeine"

    # test bad ga inputs
    test_dict["num_slots"] = -1
    assert_raises(AssertionError, test.update_inputs, test_dict)
    test_dict["num_slots"] = 100

    test_dict["init_popsize"] = 150
    assert_raises(AssertionError, test.update_inputs, test_dict)
    test_dict["init_popsize"] = 20

    test_dict["num_mevs"] = -300
    assert_raises(AssertionError, test.update_inputs, test_dict)
    test_dict["num_mevs"] = 1000

    test_dict["num_swaps"] = 20
    assert_raises(AssertionError, test.update_inputs, test_dict)
    test_dict["num_swaps"] = 1

    test_dict["num_muts"] = 15
    test.update_inputs(test_dict)
    test_dict["num_muts"] = 16
    print(test_dict)
    assert_raises(AssertionError, test.update_inputs, test_dict)
    test_dict["num_muts"] = 3

    test_dict["num_geoms"] = -2
    assert_raises(AssertionError, test.update_inputs, test_dict)
    test_dict["num_geoms"] = 3

    # num_atoms is no longer an input
    test_dict["num_atoms"] = 2
    assert_raises(InputError, test.update_inputs, test_dict)
    del test_dict["num_atoms"]

    test_dict["fit_form"] = 1
    assert_raises(AssertionError, test.update_inputs, test_dict)
    test_dict["fit_form"] = 0

    test_dict["mating_rad"] = 57
    assert_raises(AssertionError, test.update_inputs, test_dict)
    test_dict["mating_rad"] = 5

    test_dict["coef_energy"] = -5
    assert_raises(AssertionError, test.update_inputs, test_dict)
    test_dict["coef_energy"] = 0.5

    test_dict["coef_rmsd"] = -5
    assert_raises(AssertionError, test.update_inputs, test_dict)
    test_dict["coef_rmsd"] = 0.5

    # t_size no longer valid input
    test_dict["t_size"] = 25
    assert_raises(InputError, test.update_inputs, test_dict)
    del test_dict["t_size"]

    # test inchi string with slashes
    test_dict["struct_input"] = "InChI=1S/C11H22/c1-4-11(2,3)10-8-6-5-7-9-10/h10H,4-9H2,1-3H3"
    test_dict["struct_type"] = "inchi"
    test.update_inputs(test_dict)

    # cleanup
    shutil.rmtree(out_dir, ignore_errors=True)


def test_read_input():
    out_dir = os.path.join(TEST_DIR, "test_read_input")

    # remove the test directory if it already exists
    if os.path.isdir(out_dir):
        shutil.rmtree(out_dir, ignore_errors=True)
    os.mkdir(out_dir)

    # generate some output
    input_dir = {
        "struct_input": "propane",
        "output_dir": out_dir,
        "num_mevs": 5,
        "prog": "openbabel",
        "num_slots": 10,
        "init_popsize": 5,
    }
    run_kaplan(input_dir)
    # now generate a new run using the same inputs
    # make one change and start a new job
    inputs_pickle = read_input(
        os.path.join(out_dir, "kaplan_output/job_000000_propane/inputs.pickle"),
        False
    )
    # make trivial change
    inputs_pickle.num_mevs = 10
    run_kaplan(inputs_pickle)

    # unpickle the new output dir
    job1_inputs = os.path.join(out_dir, "kaplan_output/job_000001_propane/inputs.pickle")
    assert os.path.isfile(job1_inputs)
    with open(job1_inputs, "rb") as f:
        inputs1 = pickle.load(f)

    # make sure inputs are the same as those read from inputs.pickle
    for attr in inputs_pickle.__dict__:
        if attr == "num_mevs":
            # make sure change was made
            assert getattr(inputs_pickle, attr) == 10
            assert getattr(inputs1, attr) == 10
        value2 = getattr(inputs_pickle, attr)
        value3 = getattr(inputs1, attr)
        try:
            assert value2 == value3
        except ValueError:
            assert (value2 == value3).all()

    # test get_latest_job
    result = get_latest_job(out_dir)
    assert len(result) == 2
    assert result[0] == os.path.join(out_dir, "kaplan_output/job_000001_propane")
    assert result[1] == 1

    # remove the test directory as cleanup
    shutil.rmtree(out_dir, ignore_errors=True)
