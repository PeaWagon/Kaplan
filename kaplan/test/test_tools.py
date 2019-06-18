"""Test the tools module of kaplan."""

import os

from kaplan.tools import TEST_DIR, profile_function, analyse_profile


def test_profile_function():
    """Test the profile_function function from kaplan.tools.

    Notes
    -----
    This test shows the implementation for profiling a
    function. Just have to give:
    (1) Name of function to profile.
    (2) Name of dump file.
    (3) arguments followed by keyword arguments
    Use analyse_profile to convert the pstats object
    to a plain text file.

    """
    def my_adder(a, b):
        return a + b
    
    def my_list_adder(lista, listb):
        total = 0
        for i,j in zip(lista, listb):
            total += my_adder(i,j)
        return total
    
    profile_function(my_list_adder, "test_pf.dmp", range(10_000), range(10_000))
    analyse_profile("test_pf.dmp", "test_pf.log")
    # clean the profile test
    os.remove("test_pf.dmp")
    os.remove("test_pf.log")

test_profile_function()