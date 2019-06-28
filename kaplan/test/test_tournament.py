"""Test tournament module of Kaplan."""

import os

from random import seed
from copy import deepcopy

from numpy.testing import assert_raises

from kaplan.tournament import run_tournament
# , select_pmems, select_parents
from kaplan.ring import Ring, RingEmptyError
from kaplan.inputs import Inputs
from kaplan.tools import TEST_DIR


def test_run_tournament():
    """Test run_tournament function from tournament module."""
    # make tests reproducible
    seed(1234567)
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
    })
    # pre-determined tests via output of randomization
    ring = Ring(20, 5)
    run_tournament(ring, 1)
    mutations = 0
    for i, geom in enumerate(ring[3].dihedrals):
        for j, dihedral in enumerate(geom):
            if dihedral != ring[5].dihedrals[i][j]:
                mutations += 1
    assert mutations == 4
    run_tournament(ring, 2)
    mutations = 0
    for i, geom in enumerate(ring[0].dihedrals):
        for j, dihedral in enumerate(geom):
            if dihedral != ring[18].dihedrals[i][j]:
                mutations += 1
    assert mutations == 1
    mutations = 0
    for i, geom in enumerate(ring[1].dihedrals):
        for j, dihedral in enumerate(geom):
            if dihedral != ring[19].dihedrals[i][j]:
                mutations += 1
    assert mutations == 6


def test_select_pmems():
    """Test select_pmems function from tournament module."""


def test_select_parents():
    """Test select_parents function from tournament module."""
