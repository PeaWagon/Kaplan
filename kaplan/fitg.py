"""This module calculates the fitness of a
population member (pmem). The fitness is
broken into two main parts: energy and
rmsd. The energy is calculated using a
user-specified quantum chemical method
and basis set for each geometry in the pmem.
The rmsd is calculated as all the possible
pairs of rmsd between geometries."""

from kaplan.inputs import Inputs, InputError


def calc_fitness(pmem):
    """Calculate the fitness of a pmem.

    Parameters
    ----------
    pmem : object
        Pmem object containing energies,
        rmsd values, and geometries.

    Notes
    -----
    fit_form : int
        Represents the fitness formula to use.
        The only value currently available is 0,
        where fitness = CE*SE + Crmsd*Srmsd.

    Raises
    ------
    InputError
        The fit_form value is not available.

    Returns
    -------
    fitness : float

    """
    inputs = Inputs()
    if inputs.fit_form == 0:
        return abs(sum(pmem.energies))*inputs.coef_energy + \
               sum([x[2] for x in pmem.rmsds])*inputs.coef_rmsd
    raise InputError("Unsupported fitness formula.")
