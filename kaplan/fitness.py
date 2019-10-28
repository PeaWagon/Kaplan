"""

This module will be responsible for calculating
the fitness of a given pmem.

Uses
----
* pmem
* ring
* inputs

Used By
-------
* control
* tournament

IMPORTANT: Any pmem sent to this module should
already have its energy and RMSD values set,
or its fitness will be None.

There are two ways to calculate the fitness
of a pmem: with and without normalisation.

Normalisation is achieved by applying a
z-score to each energy and rmsd available
in the ring. Normalising these values
helps to make the coefficient inputs more
intuitive, as energies and rmsd values are
on a similar scale.

z-score is defined as:
z_score = (x-sample_mean)/sample_stdev

The fitness of a pmem will change over the
course of evolution (the sample mean and
sample standard deviation will change), but
its energy and rmsd values should remain
constant. Therefore, the fitness of a pmem
must be re-evaluated prior to:
(1) tournament selection,
(2) an extinction event (like plague
or deluge) that orders pmems by fitness, and
(3) the output-generation phase (when the
best pmem is chosen).

The fitness is broken into two main parts:
energy and rmsd. The energy is calculated
using a user-specified quantum chemical method
and basis set (for psi4) or using a forcefield
(for openbabel) for each geometry in the pmem.
The rmsd is calculated as all the possible
pairs of rmsd between geometries.

fit_form : int
    Represents the fitness formula to use.
    The only value currently available is 0,
    where fitness = CE*SE + Crmsd*Srmsd.


"""

from copy import deepcopy
from numpy import allclose

from kaplan.inputs import Inputs, InputError
from statistics import mean, stdev, StatisticsError


def update_all_fitness(ring):
    """A way to update all pmems from one ring.

    Notes
    -----
    If normalise is on, then updating a pmem
    one-by-one is not efficient (since the mean
    and stdev are calculated for each slot). This
    function calculates the mean and stdev once.

    """
    inputs = Inputs()
    if not inputs.normalise:
        for pmem in ring:
            if pmem:
                set_absolute_fitness(
                    pmem, inputs.fit_form,
                    inputs.coef_energy, inputs.coef_rmsd
                )
        return None

    mean_energy = ring.mean_energy
    stdev_energy = ring.stdev_energy
    mean_rmsd = ring.mean_rmsd
    stdev_rmsd = ring.stdev_rmsd

    for pmem in ring:
        if pmem:
            set_normalised_fitness(
                pmem, inputs.fit_form, inputs.coef_energy, inputs.coef_rmsd,
                mean_energy, mean_rmsd, stdev_energy, stdev_rmsd
            )


def new_pmem_fitness(ring, pmem, update_slot):
    """Calculate the fitness of a newly created pmem.

    Parameters
    ----------
    ring : kaplan.ring.Ring object
        The ring where the update is to take place.
    pmem : kaplan.pmem.Pmem object
        The new pmem.
    update_slot : int
        The location in the ring where the new pmem
        could be placed.

    Notes
    -----
    The purpose of this function is to allow calculation
    of a new pmem's fitness without it actually existing
    in the ring. This function will only really matter
    in the case where normalisation is being used, as
    the fitness of the new pmem might contribute negatively
    or positively to the mean and stdev fitness in the ring.

    Returns
    -------
    fitness : float
        The fitness of pmem.

    """
    inputs = Inputs()
    if not inputs.normalise:
        new_fitness = set_absolute_fitness(
            pmem, inputs.fit_form, inputs.coef_energy, inputs.coef_rmsd
        )
        return new_fitness

    # now determine the fitness of the pmem as if it
    # was in the ring at slot update_slot
    all_energies = ring.energies
    new_energies = [e for e in pmem.energies if e is not None]
    if ring[update_slot] is not None:
        old_energies = [e for e in ring[update_slot].energies if e is not None]
        list_replace(all_energies, old_energies, new_energies)
    else:
        all_energies += new_energies

    try:
        mean_energy = mean(all_energies)
    # ring has no valid energies
    except StatisticsError:
        mean_energy = None
    try:
        stdev_energy = stdev(all_energies)
    # only one valid energy in the ring
    except StatisticsError:
        stdev_energy = None

    # there will be no RMSDs if the number of geometries is
    # one, since no pairwise calculations can be made
    if pmem.num_geoms != 1:
        all_rmsds = ring.rmsds
        new_rmsds = [r[2] for r in pmem.rmsds if r[2] is not None]
        if ring[update_slot] is not None:
            old_rmsds = [r[2] for r in ring[update_slot].rmsds if r[2] is not None]
            list_replace(all_rmsds, old_rmsds, new_rmsds)
        else:
            all_rmsds += new_rmsds

        try:
            mean_rmsd = mean(all_rmsds)
        # ring has no valid RMSDs
        except StatisticsError:
            mean_rmsd = None
        try:
            stdev_rmsd = stdev(all_rmsds)
        # only one valid RMSD in the ring
        except StatisticsError:
            stdev_rmsd = None
    else:
        mean_rmsd = None
        stdev_rmsd = None

    fitness = set_normalised_fitness(
        pmem, inputs.fit_form, inputs.coef_energy, inputs.coef_rmsd,
        mean_energy, mean_rmsd, stdev_energy, stdev_rmsd
    )
    return fitness


def list_replace(full_list, old_values, new_values):
    """Replace values from old with new in a list.

    Parameters
    ----------
    full_list : list
        A list containing old_values.
    old_values : list
        A list of values to remove from full_list.
    new_values : list
        A list of values to be added to full_list.
        These values may already exist in the full_list,
        in which case they will be added again.

    Notes
    -----
    Can account for cases where old_values and new_values
    are not the same length. Old values MUST be in
    the original full list. Old values MUST have the same
    order in full_list as its own order. This function
    edits full_list in place.

    Returns
    -------
    full_list : list
        List after modifications.

    """
    check_none = 0
    # now update the all_values, removing those from old and adding
    # those from new
    for i, value in enumerate(full_list):
        if old_values == []:
            break
        if allclose(value, old_values[0]):
            if new_values != []:
                full_list[i] = new_values[0]
                new_values = new_values[1:]
            else:
                full_list[i] = None
                check_none += 1
            old_values = old_values[1:]
    # check if old_values has more values than new_values
    for i in range(check_none):
        full_list.remove(None)
    full_list += new_values
    return full_list


def set_fitness(ring, slot):
    """Set the fitness for a given pmem.

    Parameters
    ----------
    ring : kaplan.ring.Ring object
        The ring where the pmem lives.
    slot : int
        The slot number to update with
        a fitness value.

    Notes
    -----
    Mostly this function should be used
    when normalise is on and a new pmem
    is being considered for addition to
    the ring. Then, this function can
    calculate the fitness of the existing
    slot.

    Returns
    -------
    fitness : float
        The fitness of the pmem at the
        given slot in the ring. None can
        also be returned in the case where
        the slot is empty or the slot
        has no valid RMSD / energy values.
        None is also returned if normalise
        is being used, and there is only
        1 pmem in the ring.

    """
    pmem = ring[slot]
    if pmem is None:
        # attempted to get fitness of empty slot
        return None
    inputs = Inputs()
    if not inputs.normalise:
        fitness = set_absolute_fitness(
            pmem, inputs.fit_form, inputs.coef_energy, inputs.coef_rmsd
        )
        return fitness

    mean_energy = ring.mean_energy
    stdev_energy = ring.stdev_energy
    mean_rmsd = ring.mean_rmsd
    stdev_rmsd = ring.stdev_rmsd

    fitness = set_normalised_fitness(
        pmem, inputs.fit_form, inputs.coef_energy, inputs.coef_rmsd,
        mean_energy, mean_rmsd, stdev_energy, stdev_rmsd
    )
    return fitness


def set_normalised_fitness(pmem, fit_form, coef_energy, coef_rmsd,
                           mean_energy, mean_rmsd, stdev_energy, stdev_rmsd):
    # if all components for the fitness are None, return None
    if all(e is None for e in pmem.energies) and all(r[2] is None for r in pmem.rmsds):
        print("Warning. Pmem has no valid energies/RMSDs. Fitness is None.")
        pmem.fitness = None
        return None

    # if standard deviation is None or 0.0 (which can happen if
    # 2+ identical values are all that is in the ring)
    # then set stdev to 1 so it doesn't cause a division by
    # zero error
    if stdev_energy is None or allclose(stdev_energy, 0):
        stdev_energy = 1
    if stdev_rmsd is None or allclose(stdev_rmsd, 0):
        stdev_rmsd = 1

    # mean_energy is multiplied by -1
    # standard deviation does not change
    # since we are maximising fitness, each energy value
    # should be negated (to favour low negative energies
    # and disfavour high positive energies - in the case
    # of forcefield calculations; quantum energies will
    # always be at least zero)
    mean_energy *= -1

    # penalise bad geometries (cases where energy and RMSD are
    # None) by setting their z-score for energies and RMSDs to -4*stdev
    norm_energies = [
        (-e - mean_energy) / stdev_energy if e is not None
        else -4 * stdev_energy for e in pmem.energies
    ]
    norm_rmsds = [
        (r[2] - mean_rmsd) / stdev_rmsd if r[2] is not None
        else -4 * stdev_rmsd for r in pmem.rmsds
    ]

    # now consider that the number of rmsd values and the number of energy
    # values are not the same
    energy_sum = sum(norm_energies) / pmem.num_geoms
    if pmem.num_geoms != 1:
        rmsd_sum = sum(norm_rmsds) / pmem.num_pairs
    else:
        rmsd_sum = 0

    if fit_form == 0:
        fitness = coef_energy * energy_sum + coef_rmsd * rmsd_sum
        pmem.fitness = fitness
        return fitness
    raise InputError("No such fitness formula.")


def set_absolute_fitness(pmem, fit_form, coef_energy, coef_rmsd):
    """Calculate the absolute fitness of a pmem.

    Parameters
    ----------
    fit_form : int
        Fitness formula to use (from inputs).
    coef_energy : float
        Energy term coefficient (from inputs).
    coef_rmsd : float
        Root-mean-square deviation coefficient
        (from inputs).

    Notes
    -----
    This method calculates the fitness of a
    population member (pmem). The fitness is
    broken into two main parts: energy and
    rmsd. The energy is calculated using a
    user-specified quantum chemical method
    and basis set for each geometry in the pmem.
    The rmsd is calculated as all the possible
    pairs of rmsd between geometries.

    This function should not be called if normalise
    is being used.

    fit_form : int
        Represents the fitness formula to use.
        The only value currently available is 0,
        where fitness = CE*SE + Crmsd*Srmsd.

    Returns
    -------
    fitness : float

    """
    # if all components for the fitness are None, return None
    if all(e is None for e in pmem.energies) and all(r[2] is None for r in pmem.rmsds):
        print("Warning. Pmem has no valid energies/RMSDs. Fitness is None.")
        pmem.fitness = None
        return None
    # make sure energies is not empty
    # rmsds can be empty in the case where there
    # is only one geometry per pmem (sum of empty list
    # is zero)
    valid_energies = [e for e in pmem.energies if e is not None]
    valid_rmsds = [rmsd[2] for rmsd in pmem.rmsds if rmsd[2] is not None]

    # if any of the energies are None, then penalise the fitness according
    # to the pmem's worst value (otherwise removing the None can actually
    # improve the fitness for forcefield results since this essentially
    # sets the energy to zero - which may be lower than the real energies)
    try:
        max_energy = max(valid_energies)
    except ValueError:
        # empty sequence
        max_energy = 100
    if max_energy > 0:
        energy_sum = 0
        for energy in pmem.energies:
            if energy is not None:
                energy_sum += energy
            else:
                energy_sum += 4 * max_energy
    else:
        energy_sum = sum(valid_energies)

    if fit_form == 0:
        fitness = -1 * coef_energy * energy_sum + coef_rmsd * sum(valid_rmsds)
        pmem.fitness = fitness
        return fitness
    raise InputError("No such fitness formula.")
