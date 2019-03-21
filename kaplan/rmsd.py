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
    coords1 : np.array(shape=(n,3), dtype=float)
        The coordinates for the first
        molecule (conformer).
    coords2 : np.array(shape=(n,3), dtype=float)
        The coordinates for the second
        molecule (conformer).

    Notes
    -----
    The coordinates are given as follows:
    x coord, y coord, z coord
    Example for carbon monoxide:
    [[0.0, 0.0, 0.0],[1.0, 0.0, 0.0]]
    The distances are in atomic units.

    Returns
    -------
    rmsd : float
        The rmsd between two molecular
        geometries. This rmsd is calculated
        after applying centering and a
        rotation matrix (calculated via
        the Kabsch algorithm).

    """
    # first center each molecule
    coords1 -= rmsd.centroid(coords1)
    coords2 -= rmsd.centroid(coords2)
    # calculate the rotation matrix
    rot_matrix = rmsd.kabsch(coords1, coords2)
    # apply the rotation matrix
    coords1 = np.dot(coords1, rot_matrix)
    # finally get the rmsd
    return rmsd.rmsd(coords1, coords2)
