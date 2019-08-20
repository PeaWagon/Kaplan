"""Test the rsmd module of Kaplan."""

import os
import numpy as np

from copy import deepcopy

from vetee.job import Job, JobError
from numpy.testing import assert_raises

from kaplan.rmsd import calc_rmsd, apply_centroid
from kaplan.tools import TEST_DIR


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
    apply_centroid(mol1_coords)
    apply_centroid(mol3_coords)
    assert calc_rmsd(mol1_coords, mol3_coords) == 1.0
    # test translated/rotated hydrogen
    mol1tr = Job("xyz", os.path.join(TEST_DIR, "H2-1A-transrot.xyz"), "gaussian")
    try:
        mol1tr.setup_from_xyz()
    except JobError:
        pass
    mol1tr_coords = mol1tr.xyz_coords
    apply_centroid(mol1tr_coords)
    assert calc_rmsd(mol1_coords, mol1tr_coords) == 0.0
    # test same molecule twice
    assert calc_rmsd(mol1_coords, mol1_coords) == 0.0

    # now test ability of function to calculate the RMSD with hydrogens
    # removed
    pentanol = Job("xyz", os.path.join(TEST_DIR, "2-pentanol.xyz"), "gaussian")
    try:
        pentanol.setup_from_xyz()
    except JobError:
        pass

    pentanol_atomic_nums = [8, 6, 6, 6, 6, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

    # make sure if all atoms are removed there is a raised error
    with assert_raises(AssertionError):
        calc_rmsd(pentanol.xyz_coords, pentanol.xyz_coords,
                  atomic_nums=pentanol_atomic_nums, exclude=[8, 6, 1])

    result = calc_rmsd(pentanol.xyz_coords, pentanol.xyz_coords)
    assert np.allclose(result, 0.0)
    result = calc_rmsd(pentanol.xyz_coords, pentanol.xyz_coords,
                       atomic_nums=pentanol_atomic_nums, exclude=[1])
    assert np.allclose(result, 0.0)

    # change a hydrogen atom coordinate and then make sure carbon/oxygen rmsd
    # is still zero
    pentanol2 = deepcopy(pentanol)
    pentanol2._coords[-1][1] -= 0.5
    pentanol2._coords[-1][2] -= 0.25
    pentanol2._coords[-1][3] += 0.15

    result = calc_rmsd(pentanol2.xyz_coords, pentanol.xyz_coords)
    assert not np.allclose(result, 0.0)

    result = calc_rmsd(pentanol.xyz_coords, pentanol.xyz_coords,
                       atomic_nums=pentanol_atomic_nums, exclude=[1])
    assert np.allclose(result, 0.0)
