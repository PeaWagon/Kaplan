"""
These tests are for the ga_input module. The ga_input
module has two functions: read_ga_input and verify_ga_input.
If these tests do not pass, then it is possible that
the program's constraints have been altered.

This program uses some testfiles located in the
testfiles folder of kaplan. Specifically, in
"../Kaplan/kaplan/test/testfiles", we have:
"example_ga_input_file.txt" (a good input file)
"example2_ga_input_file.txt" (a good input file with
newlines between the inputs)
"bad1_ga_input_file.txt" (too many arguments)
"bad2_ga_input_file.txt" (too few arguments)
"bad3_ga_input_file.txt" (incorrect spelling for argument)

"""

import os

from numpy.testing import assert_raises

from kaplan.ga_input import read_ga_input, verify_ga_input


# directory for this test file
test_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'testfiles')

def test_read_ga_input():
    # good input
    read_ga_input(os.path.join(test_dir, "example_ga_input_file.txt"))
    # good input with extra spaces
    read_ga_input(os.path.join(test_dir, "example2_ga_input_file.txt"))
    # no such file error
    assert_raises(FileNotFoundError, read_ga_input, 'no-such-file')
    # num mevs appears twice
    assert_raises(ValueError, read_ga_input, os.path.join(test_dir, "bad1_ga_input_file.txt"))
    # missing num atoms
    assert_raises(ValueError, read_ga_input, os.path.join(test_dir, "bad2_ga_input_file.txt"))

def test_verify_ga_input():
    # good input with extra spaces
    ga_input_dict = read_ga_input(os.path.join(test_dir, "example2_ga_input_file.txt"))
    verify_ga_input(ga_input_dict)

    # check the good file
    #{'num_mevs': '1000', 'num_slots': '100', 'num_filled': '20',
    # 'num_geoms': '3', 'num_atoms': '10', 't_size': '7', 'num_muts': '3',
    # 'num_swaps': '1', 'pmem_dist': '5', 'fit_form': '0', 'coef_energy': '0.5',
    # 'coef_rmsd': '0.5'}
    ga_input_dict = read_ga_input(os.path.join(test_dir, "example_ga_input_file.txt"))
    verify_ga_input(ga_input_dict)

    # num geoms is spelt wrong
    ga_input_dict2 = read_ga_input(os.path.join(test_dir, "bad3_ga_input_file.txt"))
    assert_raises(ValueError, verify_ga_input, ga_input_dict2)

    ga_input_dict["num_slots"] = -1
    assert_raises(AssertionError, verify_ga_input, ga_input_dict)
    ga_input_dict["num_slots"] = 100

    ga_input_dict["num_filled"] = 150
    assert_raises(AssertionError, verify_ga_input, ga_input_dict)
    ga_input_dict["num_filled"] = 20

    ga_input_dict["num_mevs"] = -300
    assert_raises(AssertionError, verify_ga_input, ga_input_dict)
    ga_input_dict["num_mevs"] = 1000

    ga_input_dict["num_swaps"] = 20
    assert_raises(AssertionError, verify_ga_input, ga_input_dict)
    ga_input_dict["num_swaps"] = 1

    ga_input_dict["num_muts"] = 30
    assert_raises(AssertionError, verify_ga_input, ga_input_dict)
    ga_input_dict["num_muts"] = 3

    ga_input_dict["num_geoms"] = -2
    assert_raises(AssertionError, verify_ga_input, ga_input_dict)
    ga_input_dict["num_geoms"] = 3

    ga_input_dict["num_atoms"] = 2
    assert_raises(AssertionError, verify_ga_input, ga_input_dict)
    ga_input_dict["num_atoms"] = 10

    ga_input_dict["fit_form"] = 1
    assert_raises(AssertionError, verify_ga_input, ga_input_dict)
    ga_input_dict["fit_form"] = 0

    ga_input_dict["pmem_dist"] = 57
    assert_raises(AssertionError, verify_ga_input, ga_input_dict)
    ga_input_dict["pmem_dist"] = 5

    ga_input_dict["coef_energy"] = -5
    assert_raises(AssertionError, verify_ga_input, ga_input_dict)
    ga_input_dict["coef_energy"] = 0.5

    ga_input_dict["coef_rmsd"] = -5
    assert_raises(AssertionError, verify_ga_input, ga_input_dict)
    ga_input_dict["coef_rmsd"] = 0.5

    ga_input_dict["t_size"] = 25
    assert_raises(AssertionError, verify_ga_input, ga_input_dict)

if __name__ == "__main__":
    test_read_ga_input()
    test_verify_ga_input()
