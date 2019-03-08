"""This module calculates the fitness of a
population member (pmem). The fitness is
broken into two main parts: energy and
rmsd. The energy is calculated using a
user-specified quantum chemical method
and basis set for each geometry in the pmem.
The rmsd is calculated as all the possible
pairs of rmsd between geometries."""

from math import factorial

import numpy as np

from kaplan.energy import run_energy_calc
from kaplan.rmsd import calc_rmsd
from kaplan.inputs import Inputs

# TODO incorporate parser attribute "prog"
# (program) such that a user could specify
# another program (other than psi4) to
# calculate energies


def sum_energies(xyz_coords):
    """Sum the energy calculations for a pmem.

    Parameters
    ----------
    xyz_coords : list
        The length of xyz_coords is how many
        geometries are in the pmem. Each geometry
        is specified as:
        np.array([[x1,y1,z1], [x2,y2,z2],
                    ..., [xn,yn,zn]])
        The x,y,z coordinates for each atom
        (for a total of n atoms). The coordinates
        are given as floating point numbers in
        atomic units.

    """
    energies = np.zeros(len(xyz_coords), float)
    for i, xyz in enumerate(xyz_coords):
        try:
            energies[i] = run_energy_calc(xyz)
        except Exception:
            # if there is a convergence error (atom too close)
            # give an energy of zero
            print("Warning: non-convergence for molecule.")
            energies[i] = 0
    return abs(sum(energies))


def sum_rmsds(xyz_coords):
    """Sum the rmsd calculations for a pmem.

    Parameters
    ----------
    xyz_coords : list
        The length of xyz_coords is how many
        geometries are in the pmem. Each geometry
        is specified as:
        np.array([[x1,y1,z1], [x2,y2,z2],
                    ..., [xn,yn,zn]])
        The x,y,z coordinates for each atom
        (for a total of n atoms). The coordinates
        are given as floating point numbers in
        atomic units.

    """
    num_geoms = len(xyz_coords)
    # n choose k = n!/(k!(n-k)!)
    num_pairs = int(factorial(num_geoms)/(2*factorial(num_geoms-2)))
    rmsd_values = np.zeros(num_pairs, float)
    pairs = all_pairs_gen(len(xyz_coords))
    for i in range(num_pairs):
        ind1, ind2 = next(pairs)
        result = calc_rmsd(xyz_coords[ind1], xyz_coords[ind2])
        rmsd_values[i] = calc_rmsd(xyz_coords[ind1], xyz_coords[ind2])
    return sum(rmsd_values)


def all_pairs_gen(num_geoms):
    """Yield indices of two geometries.

    Note
    ----
    This is a generator function.

    """
    for i in range(num_geoms-1):
        for j in range(i+1, num_geoms):
            yield (i, j)


def calc_fitness(all_coords):
    """Calculate the fitness of a pmem.

    Parameters
    ----------
    all_coords : list
        Each member of the list is:
        np.array(shape=(n,3), dtype=float)
        Where n is the number of atoms. These
        arrays specify the xyz coordinates of
        one conformer. The length of the list
        is equal to num_geoms.

    Notes
    -----
    fit_form : int
        Represents the fitness formula to use.
        The only value currently available is 0,
        where fitness = CE*SE + Crmsd*Srmsd.

    Raises
    ------
    ValueError
        The fit_form value is not available.

    Returns
    -------
    fitness : float

    """
    inputs = Inputs()
    if inputs.fit_form == 0:
        return sum_energies(all_coords)*inputs.coef_energy + sum_rmsds(all_coords)*inputs.coef_rmsd
    raise ValueError("Unsupported fitness formula.")
