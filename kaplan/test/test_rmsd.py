"""Test the rsmd module of Kaplan."""

import os
import shutil
import numpy as np

from copy import deepcopy

from numpy.testing import assert_raises

from kaplan.rmsd import calc_rmsd, apply_centroid
from kaplan.tools import TEST_DIR
from kaplan.geometry import create_obmol, get_coords, get_atomic_nums


def test_calc_rmsd():
    """Test the calc_rmsd function from the rsmd module."""
    # test the A1 and A3 hydrogen molecules
    out_dir = os.path.join(TEST_DIR, "test_calc_rmsd")
    # remove the test directory if it already exists
    if os.path.isdir(out_dir):
        shutil.rmtree(out_dir, ignore_errors=True)
    os.mkdir(out_dir)

    test = create_obmol(os.path.join(TEST_DIR, "H2-1A.xyz"), "xyz", True)
    mol1_coords = get_coords(test)
    apply_centroid(mol1_coords)

    test = create_obmol(os.path.join(TEST_DIR, "H2-3A.xyz"), "xyz", True)
    mol3_coords = get_coords(test)
    apply_centroid(mol3_coords)

    assert calc_rmsd(mol1_coords, mol3_coords) == 1.0

    # test translated/rotated hydrogen
    test = create_obmol(os.path.join(TEST_DIR, "H2-1A-transrot.xyz"), "xyz", True)
    mol1tr_coords = get_coords(test)
    apply_centroid(mol1tr_coords)

    assert calc_rmsd(mol1_coords, mol1tr_coords) == 0.0
    # test same molecule twice
    assert calc_rmsd(mol1_coords, mol1_coords) == 0.0

    # now test ability of function to calculate the RMSD with hydrogens
    # removed
    test = create_obmol(os.path.join(TEST_DIR, "2-pentanol.xyz"), "xyz", True)

    pentanol_atomic_nums = [8, 6, 6, 6, 6, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    atomic_nums = get_atomic_nums(test)
    assert atomic_nums == pentanol_atomic_nums
    test_coords = get_coords(test)

    # make sure if all atoms are removed there is a raised error
    with assert_raises(AssertionError):
        calc_rmsd(
            test_coords, test_coords,
            atomic_nums=atomic_nums, exclude=[8, 6, 1]
        )

    result = calc_rmsd(test_coords, test_coords)
    assert np.allclose(result, 0.0)
    result = calc_rmsd(
        test_coords, test_coords,
        atomic_nums=pentanol_atomic_nums, exclude=[1]
    )
    assert np.allclose(result, 0.0)

    # change a hydrogen atom coordinate and then make sure carbon/oxygen rmsd
    # is still zero
    pentanol2 = deepcopy(test_coords)
    pentanol2[-1][0] -= 0.5
    pentanol2[-1][1] -= 0.25
    pentanol2[-1][2] += 0.15

    result = calc_rmsd(pentanol2, test_coords)
    assert not np.allclose(result, 0.0)

    result = calc_rmsd(
        test_coords, test_coords,
        atomic_nums=pentanol_atomic_nums, exclude=[1]
    )
    assert np.allclose(result, 0.0)

    # cleanup
    shutil.rmtree(out_dir, ignore_errors=True)
