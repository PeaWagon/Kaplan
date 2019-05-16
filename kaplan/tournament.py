"""The tournament module decides which pmems to
pick from the ring in order to apply updates
to the population."""

from math import inf

from random import choice

from kaplan.inputs import Inputs
from kaplan.ring import RingEmptyError
from kaplan.mutations import generate_children, single_parent_mutation


def run_tournament(ring, current_mev):
    """Run the tournament (i.e. a mating event).

    Parameters
    ----------
    ring : object
        Ring object.
    current_mev : int
        The current mating event number. Used
        to give pmems birthdays.

    Returns
    -------
    None

    """
    inputs = Inputs()
    # check ring has enough pmems for a tournament
    if not ring.num_filled:
        raise RingEmptyError("Need at least one pmem to run a tournament.")

    # choose a random pmem for a tournament
    # any pmem in its mating radius is considered
    selected_pmems = select_pmems(inputs.pmem_dist, ring)

    # select parents and worst slots by fitness
    parents, worst = select_parents(selected_pmems, ring)

    parent1 = ring[parents["p1"][1]].dihedrals
    try:
        parent2 = ring[parents["p2"][1]].dihedrals
    # there is only one pmem in the mating radius
    except AttributeError:
        child1 = single_parent_mutation(parent1, inputs.num_muts)
        child2 = None
    else:
        # generate children
        child1, child2 = generate_children(
                             parent1, parent2,
                             inputs.num_muts, inputs.num_swaps,
                             inputs.num_cross
                         )

    # put children in ring
    ring.update(child1, worst["w1"][1], current_mev)
    if child2:
        ring.update(child2, worst["w2"][1], current_mev)


def select_pmems(pmem_dist, ring):
    """Randomly selected pmems.

    Parameters
    ----------
    pmem_dist : int
        How many slots (to the left and right) to
        consider when selecting pmems from the ring.
    ring : object
        The ring from which to pick the pmems.

    """
    # get a list of indices representing filled slots
    occupied = [i for i in range(ring.num_slots) if ring[i] is not None]
    # from the occupied slots, choose a random pmem as parent1
    parent1 = choice(occupied)
    selection = ring.mating_radius(parent1, pmem_dist)
    return selection


def select_parents(selected_pmems, ring):
    """Sort pmems by fitness values, choose best 2.

    Parameters
    ----------
    selected_pmems : list
        The indices of the ring that are in the tournament.

    ring : object
        Instance of Ring class.

    Returns
    -------
    tuple : int, int
        Indices of two best pmems to be used as parents.

    """
    # values are tuple(fitness, ring index)
    best_two = {"p1": (-inf, None), "p2": (-inf, None)}
    worst_two = {"w1": (inf, None), "w2": (inf, None)}
    for s in select_pmems:
        try:
            fit = ring[s].fitness
        # ring at slot is empty
        # AttributeError: 'NoneType' object has no attribute 'fitness'
        except AttributeError:
            fit = -inf
        if fit > best_two["p1"][0]:
            best_two["p1"] = (fit, s)
        elif fit > best_two["p2"][0]:
            best_two["p2"] = (fit, s)
        elif fit < worst_two["w1"][0]:
            worst_two["p2"] = (fit, s)
        elif fit < worst_two["w2"][0]:
            worst_two["p2"] = (fit, s)
    return best_two, worst_two
