"""This module is repsonsible for calculating the rmsd
(root-mean square deviation) between two sets of
coordinates. It uses the rmsd library, which can
be found here: https://github.com/charnley/rmsd."""

import numpy as np
import rmsd


def calc_rmsd(coords1, coords2):
    """Calculate root-mean square deviation using kabsch algorithm.

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
    The distances are in angstroms.
    This function assumes the two sets
    of coordinates have already been centred
    by using the apply_centroid function.

    Returns
    -------
    rmsd : float
        The rmsd between two molecular
        geometries. This rmsd is calculated
        after applying a rotation matrix
        (calculated via the Kabsch algorithm).
        No centering is applied.

    """
    # calculate the rotation matrix
    rot_matrix = rmsd.kabsch(coords1, coords2)
    # apply the rotation matrix
    coords1 = np.dot(coords1, rot_matrix)
    # finally get the rmsd
    return rmsd.rmsd(coords1, coords2)


def apply_centroid(coords):
    """Get coordinates after applying a centroid translation.

    Parameters
    ----------
    coords : np.array((num_atoms, 3), float)
        The input coordinates as a numpy array.

    Notes
    -----
    Same as in rmsd python package. Centroid is the mean position:
    centroid_x = sum(coords_x)/len(coords_x) (where it is the same
    for y and z coordinates).

    Returns
    -------
    np.array((num_atoms,3), float) centred about centroid coordinate.

    """
    coords -= coords.mean(axis=0)
    return coords
