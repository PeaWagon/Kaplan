

import os

from pprint import pprint
from kaplan.ga_input import read_ga_input, verify_ga_input

from numpy.testing import assert_raises

# directory for this test file
test_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'testfiles')

def test_read_ga_input():
    ga_input_dict = read_ga_input(os.path.join(test_dir, "example_ga_input_file.txt"))
    assert_raises(FileNotFoundError, read_ga_input, 'no-such-file')


def test_verify_ga_input():
    ga_input_dict = read_ga_input(os.path.join(test_dir, "example_ga_input_file.txt"))
    verify_ga_input(ga_input_dict)

if __name__ == "__main__":
    test_read_ga_input()
    test_verify_ga_input()
