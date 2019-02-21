"""Test tournament module of Kaplan."""

import os

from vetee.xyz import Xyz
from numpy.testing import assert_raises

from kaplan.tournament import run_tournament
# , select_pmems, select_parents
from kaplan.ring import Ring, RingEmptyError
from kaplan.inputs import Inputs

# directory for this test file
TEST_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'testfiles')


def test_run_tournament():
    """Test run_tournament function from tournament module."""
    # ring init:
    # num_geoms, num_atoms, num_slots, pmem_dist, fit_form,
    # coef_energy, coef_rmsd, parser
    inputs = Inputs()
    inputs.update_inputs({
        "struct_input": os.path.join(TEST_DIR, "1,3-butadiene.xyz"),
        "struct_type": "xyz",
        "charge": 0,
        "multip": 1,
        "num_geoms": 3,
        "num_slots": 20,
        "pmem_dist": 2
    })
    ring = Ring()
    ring.fill(3, 0)
    run_tournament(ring, 0)


def test_select_pmems():
    """Test select_pmems function from tournament module."""


def test_select_parents():
    """Test select_parents function from tournament module."""
