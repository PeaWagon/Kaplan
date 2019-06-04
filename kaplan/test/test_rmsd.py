"""Test the rsmd module of Kaplan."""

import os
import numpy as np

from vetee.job import Job, JobError

from kaplan.rmsd import calc_rmsd


# directory for this test file
TEST_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'testfiles')


def test_calc_rmsd():
    """Test the calc_rmsd function from the rsmd module."""
    # test the A1 and A3 hydrogen molecules
    mol1 = Job("xyz", os.path.join(TEST_DIR, "H2-1A.xyz"), "gaussian")
    try:
        mol1.setup_from_xyz()
    except JobError:
        pass
    mol3 = Job("xyz", os.path.join(TEST_DIR, "H2-3A.xyz"), "gaussian")
    try:
        mol3.setup_from_xyz()
    except JobError:
        pass
    
    mol1_coords = mol1.xyz_coords
    mol3_coords = mol3.xyz_coords
    assert calc_rmsd(mol1_coords, mol3_coords) == 1.0
    # test translated/rotated hydrogen
    mol1tr = Job("xyz", os.path.join(TEST_DIR, "H2-1A-transrot.xyz"), "gaussian")
    try:
        mol1tr.setup_from_xyz()
    except JobError:
        pass
    mol1tr_coords = mol1tr.xyz_coords
    assert calc_rmsd(mol1_coords, mol1tr_coords) == 0.0
    # test same molecule twice
    assert calc_rmsd(mol1_coords, mol1_coords) == 0.0
