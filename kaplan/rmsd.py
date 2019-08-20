"""This module is repsonsible for calculating the rmsd
(root-mean square deviation) between two sets of
coordinates. It uses the rmsd library, which can
be found here: https://github.com/charnley/rmsd."""

import numpy as np
import rmsd


def calc_rmsd(coords1, coords2, atomic_nums=None, exclude=None):
    """Calculate root-mean square deviation using kabsch algorithm.

    Parameters
    ----------
    coords1 : np.array(shape=(n,3), dtype=float)
        The coordinates for the first
        molecule (conformer).
    coords2 : np.array(shape=(n,3), dtype=float)
        The coordinates for the second
        molecule (conformer).
    atomic_nums : list(int)
        The atomic numbers for each atom in
        the molecule for coords1 and coords2.
        The atom ordering for both sets of
        coordinates should be the same. Defaults
        to None. Is only required if exlude is not
        None.
    exclude : list(int)
        List of atomic numbers to exclude from
        the RMSD calculation. Defaults to None.
        If None, all atomic numbers are included.

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
    if exclude is None:
        # calculate the rotation matrix
        rot_matrix = rmsd.kabsch(coords1, coords2)
        # apply the rotation matrix
        coords1 = np.dot(coords1, rot_matrix)
        # finally get the rmsd
        return rmsd.rmsd(coords1, coords2)

    # determine how many atoms to exclude
    assert atomic_nums is not None
    # make sure not to doubly count any atomic numbers
    exclude = set(exclude)
    orig_atom_count = len(coords1)
    new_atom_count = orig_atom_count
    assert len(atomic_nums) == orig_atom_count
    for exclude_atom in exclude:
        new_atom_count -= atomic_nums.count(exclude_atom)
    # make sure the RMSD is calculated for at least
    # one atom type
    assert new_atom_count != 0

    # create new coordinates for the geometry
    # without the excluded atoms
    new_coords1 = np.zeros((new_atom_count, 3), float)
    new_coords2 = np.zeros((new_atom_count, 3), float)

    new_coord_index = 0
    for i, atom in enumerate(coords1):
        if atomic_nums[i] not in exclude:
            new_coords1[new_coord_index] = atom
            new_coords2[new_coord_index] = coords2[i]
            new_coord_index += 1

    # calculate the rotation matrix
    rot_matrix = rmsd.kabsch(new_coords1, new_coords2)
    # apply the rotation matrix
    new_coords1 = np.dot(new_coords1, rot_matrix)
    # finally get the rmsd
    return rmsd.rmsd(new_coords1, new_coords2)


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
