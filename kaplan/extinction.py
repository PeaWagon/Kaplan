"""Kaplan extinction module.

This module provides functions that operate on ring
objects. These optional functions are designed to
reset the ring to a more exploratory state.

"""

from kaplan.ring import RingEmptyError
from kaplan.inputs import InputError
from kaplan.tournament import sortby

from random import uniform, randint

def apply_extinction(ring, operator, normalise):
    """Apply an extinction event to a ring object.

    Parameters
    ----------
    ring : kaplan.ring.Ring object
        Contains information about pmems.
    operator : str
        Name of extinction operator to apply.
        Should be one of: asteroid, plague, agathic,
        or deluge.
    normalise : bool
        From inputs.normalise. If True, then operators
        that are applied based on fitness (such as deluge
        and plague) will re-evaluate all fitness values
        prior to the extinction event. If False, then
        the operators use the absolute fitness value.
    
    Raises
    ------
    NotImplementedError
        The operator is not available (or was
        misspelled).
    RingEmptyError
        Ring contains no pmems and thus extinction
        will do nothing.
        Also raised if deluge/plague is called with
        a population where no fitness values have
        been initialised (cannot take max of None).
    InputError
        One or more pmems have slot numbers that
        are out of bounds for the ring.
    AssertionError
        Number of occupied slots is not the same
        as the num_filled ring attribute.

    Notes
    -----
    The extinction operators can be called directly,
    but this function will check that the inputs are
    valid first (which makes it a safer and more
    computationally expensive function).

    Returns
    -------
    ring object after applying extinction operator.

    """
    # check ring is not already empty and that occupied matches num_filled
    occupied_indices = ring.occupied
    assert len(occupied_indices) == ring.num_filled
    if occupied_indices == []:
        raise RingEmptyError(f"Ring is empty. Cannot apply {operator}.")

    # check pmems are all within allowed boundaries
    if any(slot not in range(ring.num_slots) for slot in occupied_indices):
        raise InputError("Ring contains pmems with ring_loc outside of ring.")

    # in case user inputs operator with uppercase
    operator = operator.lower()

    # apply extinction operator
    if operator == "deluge":
        deluge(ring, occupied_indices, normalise)
    elif operator == "plague":
        plague(ring, normalise)
    elif operator == "asteroid":
        asteroid(ring, occupied_indices)
    elif operator == "agathic":
        agathic(ring)
    else:
        raise NotImplementedError(f"No such extinction operator: {operator}")
    
    return ring
    

def asteroid(ring, occupied_indices):
    """Apply the asteroid extinction operator to the ring.
    
    The asteroid operator deletes a contiguous segment
    of the ring. The size of the deleted segment is chosen
    at random between 10-90% of the total number of slots.
    The starting slot of the segment (for deletion) is chosen
    at random, provided it does not kill the entire population,
    in which case another starting slot is chosen.
    
    """
    while True:
        start_slot = randint(0, ring.num_slots-1)
        asteroid_size = int(uniform(0.1, 0.9)*ring.num_slots)
        end_slot = start_slot + asteroid_size - 1
        
        #print("start slot", start_slot)
        #print("asteroid size", asteroid_size)
        #print("end slot", end_slot)
        # no ring wrapping
        if end_slot < ring.num_slots:
            #print("No wrap")
            if any(slot < start_slot or slot > end_slot for slot in occupied_indices):
                for slot in range(start_slot, end_slot+1):
                    ring[slot] = None
                return ring
        
        # consider ring wrapping
        else:
            #print("Wrap")
            overflow = end_slot - ring.num_slots
            # overflow is the last index of pmem to be killed
            #print("overflow", overflow)
            if any(slot > overflow and slot < start_slot for slot in occupied_indices):
                for slot in range(overflow+1):
                    ring[slot] = None
                for slot in range(start_slot, ring.num_slots):
                    ring[slot] = None
                return ring


def plague(ring, normalise):
    """Apply the plague extinction operator to the ring.
    
    Removes the fraction of the population with the lowest
    fitness, where the fraction is chosen as 10-90% of
    the current population size.
    
    """
    if normalise:
        for pmem in ring.occupied:
            ring.set_fitness(ring[pmem])
    # essentially the same procedure to sort pmems
    # as in the tournament module
    index_fit_pairs = sortby(ring, "fitness")
    # remove empty pmems from pairs
    empty = ring.num_slots - ring.num_filled
    index_fit_pairs = index_fit_pairs[empty:]

    # choose random amount of pmems to kill
    plague_size = int(uniform(0.1, 0.9)*ring.num_filled)
    # apply plague operator
    for death in range(plague_size):
        slot = index_fit_pairs[death][0]
        ring[slot] = None
    return ring



def agathic(ring):
    """Apply the agathic extinction operator to the ring.
    
    The agathic operator works in the same way as the
    plague operator, except the pmems are sorted by
    age and not fitness. The oldest fraction are removed.
    The fraction to remove is 10-90% of the current
    population size.
    
    """
    # essentially the same procedure to sort pmems
    # as in the tournament module
    index_bday_pairs = sortby(ring, "birthday")
    # remove empty pmems from pairs
    empty = ring.num_slots - ring.num_filled
    index_bday_pairs = index_bday_pairs[empty:]

    # choose random amount of pmems to kill
    agathic_size = int(uniform(0.1, 0.9)*ring.num_filled)
    # apply plague operator
    for death in range(agathic_size):
        slot = index_bday_pairs[death][0]
        ring[slot] = None
    return ring


def deluge(ring, occupied_indices, normalise):
    """Apply the deluge extinction operator to the ring.
    
    This operator first chooses a water-level, which is
    a fraction of the maximum (current) population fitness.
    The water-level is randomly chosen between 10-90%
    inclusive. Any pmem in the ring with a fitness less
    than the water-level is killed (deleted).

    For example, if the fraction is randomly chosen at 20%
    and the maximum fitness is 10, then any pmem with fitness
    less than 2 is removed from the population.

    """
    # first make sure existing pmems have an up-to-date fitness
    if normalise:
        for pmem in ring.occupied:
            ring.set_fitness(ring[pmem])
    # make sure maximum fitness is calculable
    _, max_fit = ring.best_pmem
    if max_fit is None:
        raise RingEmptyError("Need at least one pmem with a fitness value to call deluge.")
    # if max fitness is negative, don't want to kill all pmems
    # instead treat as a minimisation problem
    if max_fit < 0:
        water_level = uniform(1.1, 1.9)*max_fit*-1
        for slot in occupied_indices:
            if ring[slot].fitness is None or -1*ring[slot].fitness > water_level:
                ring[slot] = None
        return ring

    # in the case where the fitness values are positive, then it should be
    # possible to kill all pmems except those with at least 90% of the max fitness
    water_level = uniform(0.1, 0.9)*max_fit
    for slot in occupied_indices:
        if ring[slot].fitness is None or ring[slot].fitness < water_level:
            ring[slot] = None
    return ring
