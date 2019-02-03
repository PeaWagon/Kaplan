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

from kaplan.energy import run_energy_calc, prep_psi4_geom
from kaplan.rmsd import calc_rmsd

# TODO incorporate parser attribute "prog"
# (program) such that a user could specify
# another program (other than psi4) to
# calculate energies


def sum_energies(xyz_coords, charge, multip, method, basis):
    """Sum the energy calculations for a pmem.

    Parameters
    ----------
    xyz_coords : list
        The length of xyz_coords is how many
        geometries are in the pmem. Each geometry
        is specified as:
        [["A",x1,y1,z1], ["B",x2,y2,z2],
                    ..., ["C",xn,yn,zn]]
        Where the letters are the elements, and
        there are x,y,z coordinates for each atom
        (for a total of n atoms). The coordinates
        are given as integers.
    charge : int
        The charge of the molecule.
    multip : int
        The multiplicity of the molecule.
    method : str
        The quantum chemical method to use to
        calculate the energy.
    basis : str
        The basis set to use to calculate the
        energy.

    """
    energies = np.zeros(len(xyz_coords), float)
    for i, xyz in enumerate(xyz_coords):
        try:
            energies[i] = run_energy_calc(prep_psi4_geom(xyz, charge, multip), method, basis)
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
        [["A",x1,y1,z1], ["B",x2,y2,z2],
                    ..., ["C",xn,yn,zn]]
        Where the letters are the elements, and
        there are x,y,z coordinates for each atom
        (for a total of n atoms). The coordinates
        are given as integers.

    """
    num_geoms = len(xyz_coords)
    rmsd_values = np.zeros(len(xyz_coords), float)
    # n choose k = n!/(k!(n-k)!)
    num_pairs = int(factorial(num_geoms)/(2*factorial(num_geoms-2)))
    pairs = all_pairs_gen(len(xyz_coords))
    for i in range(num_pairs):
        ind1, ind2 = next(pairs)
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


def calc_fitness(fit_form, sum_energy, coef_energy, sum_rmsd, coef_rmsd):
    """Calculate the fitness of a pmem.

    Parameters
    ----------
    fit_form : int
        Represents the fitness formula to use.
        The only value currently available is 0,
        where fitness = CE*SE + Crmsd*Srmsd.
    sum_energy : float
        The summation of all of the individual
        energy calculations for each of the geometries.
    coef_energy : float
        The energy coefficient in the fitness formula.
    sum_rmsd : float
        The summation of all rmsd when comparing
        pairs of geometries.
    coef_rmsd : float
        The rmsd coefficient in the fitness formula.

    Raises
    ------
    ValueError
        The fit_form value is not available.

    Returns
    -------
    fitness : float

    """
    if fit_form == 0:
        return sum_energy*coef_energy + sum_rmsd*coef_rmsd
    raise ValueError("Unsupported fitness formula.")
