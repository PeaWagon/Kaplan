"""The tournament module decides which pmems to
pick from the ring in order to apply updates
to the population."""

import numpy as np

from kaplan.inputs import Inputs
from kaplan.ring import RingEmptyError
from kaplan.mutations import generate_children


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
    if inputs.t_size > ring.num_filled:
        raise RingEmptyError("Not enough pmems to run a tournament.")

    # choose random slots for a tournament
    selected_pmems = select_pmems(inputs.t_size, ring)

    # select parents by fitness
    parents = select_parents(selected_pmems, ring)

    parent1 = ring[parents[0]].dihedrals
    parent2 = ring[parents[1]].dihedrals

    # generate children
    children = generate_children(parent1, parent2, inputs.num_muts, inputs.num_swaps)

    # put children in ring
    ring.update(parents[0], children[1], current_mev)
    ring.update(parents[1], children[0], current_mev)


def select_pmems(number, ring):
    """Randomly selected pmems.

    Parameters
    ----------
    number : int
        How many pmems to pick.
    ring : object
        The ring from which to pick the pmems.

    """
    # get a list of indices representing filled slots
    occupied = [i for i in range(ring.num_slots) if ring[i] is not None]
    # from the occupied slots, choose number of them
    # without replacement (i.e. don't pick the same
    # pmem twice)
    selection = np.random.choice(occupied, number, replace=False)
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
    fit_vals = np.array([ring[i].fitness for i in selected_pmems])
    # from here:
    # https://stackoverflow.com/questions/6910641/how-do-i-get-indices-of-n-maximum-values-in-a-numpy-array
    # use numpy to get the two best fitness value indices
    # from the list and link it to ring index
    # note: this makes generator object
    parents_gen = (selected_pmems[parent] for parent in np.argpartition(fit_vals, -2)[-2:])
    parents = [next(parents_gen)]
    parents.append(next(parents_gen))
    return parents
