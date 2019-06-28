"""This test runs the main Kaplan program
with some dummy input files."""

import os

from kaplan.control import run_kaplan
from kaplan.tools import TEST_DIR


def test_run_kaplan():
    """Test run_kaplan function from control module."""
    
    test1 = {
        "struct_input": "propane",
        "struct_type": "name",
        "charge": 0,
        "multip": 1,
        "num_mevs": 20,
        "init_popsize": 5,
        "t_size": 5,
        "num_slots": 20
    }
    #run_kaplan(test1)
    
    # caffeine might be too hard of a test right now
    test2 = {
        "struct_input": os.path.join(TEST_DIR, "caffeine.xyz"),
        "struct_type": "xyz",
        "charge": 0,
        "multip": 1,
        "num_mevs": 20,
        "init_popsize": 5,
        "t_size": 5,
        "num_slots": 20
    }
    #run_kaplan(test2)

    test3 = {
        "struct_input": "threonine",
        "struct_type": "name",
        "charge": 0,
        "multip": 1,
        "num_mevs": 20,
        "init_popsize": 5,
        "t_size": 5,
        "num_slots": 20 
    }
    run_kaplan(test3)
