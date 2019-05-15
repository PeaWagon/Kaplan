"""Test geometry module from Kaplan."""

import os

from kaplan.geometry import get_geom_from_dihedrals
from kaplan.inputs import Inputs

# directory for this test file
TEST_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'testfiles')


def test_get_geom_from_dihedrals():
    """Test get_geom_from_dihedrals function of geometry module."""
    pass
