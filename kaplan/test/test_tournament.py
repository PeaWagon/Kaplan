"""Test tournament module of Kaplan."""

import os

from random import seed
from copy import deepcopy

from numpy.testing import assert_raises

from kaplan.tournament import run_tournament
from kaplan.ring import Ring, RingEmptyError
from kaplan.inputs import Inputs
from kaplan.tools import TEST_DIR


def test_run_tournament():
    """Test run_tournament function from tournament module."""
    inputs = Inputs()
    inputs.update_inputs({
        "struct_input": os.path.join(TEST_DIR, "1,3-butadiene.xyz"),
        "struct_type": "xyz",
        "charge": 0,
        "multip": 1,
        "num_geoms": 3,
        "init_popsize": 5,
        "num_slots": 20,
        "mating_rad": 2,
        "num_cross": 0,
        "num_muts": 6,
        "num_swaps": 0,
        "normalise": False,
    })

    def print_this():
        for i, slot in enumerate(ring):
            if slot is None:
                print(f"Slot {i} Empty")
            else:
                print(f"Slot {i} Fitness: {slot.fitness}")
        print("-" * 50)

    ring = Ring(20, 5)
    print_this()
    for i in range(10):
        run_tournament(ring, i, True)
        print_this()

    inputs.update_inputs({
        "struct_input": os.path.join(TEST_DIR, "1,3-butadiene.xyz"),
        "struct_type": "xyz",
        "charge": 0,
        "multip": 1,
        "num_geoms": 3,
        "init_popsize": 5,
        "num_slots": 20,
        "mating_rad": 2,
        "num_cross": 0,
        "num_muts": 6,
        "num_swaps": 0,
        "normalise": True,
    })

    ring = Ring(20, 5)
    print_this()
    for i in range(10):
        run_tournament(ring, i, True)
        print_this()
