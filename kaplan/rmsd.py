"""This module is repsonsible for calculating the rmsd
(root-mean square deviation) between two sets of
coordinates. It uses the rmsd library, which can
be found here: https://github.com/charnley/rmsd."""

import numpy as np
import rmsd


def calc_rmsd(coords1, coords2):
    """Calculate root-mean square deviation.

    Parameters
    ----------
    coords1 : list(list(str,int,int,int))
        The coordinates for the first
        molecule (conformer).
    coords2 : list(list(str,int,int,int))
        The coordinates for the second
        molecule (conformer).

    Notes
    -----
    The coordinates are given as follows:
    atom name, x coord, y coord, z coord
    Example for carbon monoxide:
    [[C, 0.0, 0.0, 0.0],[O, 1.0, 0.0, 0.0]]
    The distances are in Angstroms.

    Returns
    -------
    rmsd : float
        The rmsd between two molecular
        geometries. This rmsd is calculated
        after applying centering and a
        rotation matrix (calculated via
        the Kabsch algorithm).

    """
    # trivial check
    assert len(coords1) == len(coords2)
    # first turn lists into numpy arrays (otherwise
    # cannot perform some rmsd operations)
    # and excise atom names (not needed for rmsd)
    mol1 = np.array([[coords1[i][1], coords1[i][2], coords1[i][3]] for i in range(len(coords1))])
    mol2 = np.array([[coords2[i][1], coords2[i][2], coords2[i][3]] for i in range(len(coords2))])
    # first center each molecule
    mol1 -= rmsd.centroid(mol1)
    mol2 -= rmsd.centroid(mol2)
    # calculate the rotation matrix
    rot_matrix = rmsd.kabsch(mol1, mol2)
    # apply the rotation matrix
    mol1 = np.dot(mol1, rot_matrix)
    # finally get the rmsd
    return rmsd.rmsd(mol1, mol2)
