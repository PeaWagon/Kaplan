

import numpy as np
from random import sample, randint
from copy import deepcopy

# values for dihedral angles in degrees
MIN_VALUE = 0
MAX_VALUE = 360

"""
Might want to add another mutate function
that changes the dihedral angle by +/- a set
value (to avoid massive angular changes).

"""

def generate_children(parent1, parent2, num_muts, num_swaps):
    """Make some new pmems for the ring.

    Parameters
    ----------
    parent1 : list(list(int))
        list of lists of dihedral angles
    parent2 : list(list(int))
        list of lists of dihedral angles
    num_muts : int
        maximum number of mutations to perform
    num_swaps : int
        maximum number of swaps to perform

    Returns
    -------
    two new sets of dihedral angles with which
    to make pmem objects (list(list(int))x2

    """
    # check parent sizes
    assert len(parent1) == len(parent2)
    # check num_swaps
    assert num_swaps <= len(parent1)
    # check num_muts
    assert num_muts <= len(parent1[0])
    child1 = deepcopy(parent1)
    child2 = deepcopy(parent2)
    child1, child2 = swap(child1, child2, num_swaps)
    child1 = mutate(child1, num_muts)
    child2 = mutate(child2, num_muts)
    return child1, child2

def mutate(dihedrals, num_muts):
    """Mutate a list of dihedral angles.

    Parameters
    ----------
    dihedrals : list(list (int))
        List of lists of dihedral angles to mutate.
    num_muts : int
        Maximum number of mutations to perform on list.

    Returns
    -------
    dihedrals : list(list (int))
        Dihedral angles after mutations.

    """
    for geom in range(len(dihedrals)):
        # choose how many mutations to do
        num_muts = randint(0, num_muts)
        # choose where to do the mutations
        mut_ind = sample(range(len(dihedrals[0])), num_muts)
        for m in mut_ind:
            dihedrals[geom][m] = randint(MIN_VALUE, MAX_VALUE-1)
    return dihedrals

def swap(child1, child2, num_swaps):

    # choose how many swaps to do
    num_swaps = randint(0, num_swaps)
    # choose where to do the swaps
    swap_ind = sample(range(len(child1)), num_swaps)
    # apply swaps
    for s in swap_ind:
        child1_value = child1[s]
        child1[s] = child2[s]
        child2[s] = child1_value
    # return updated pmems
    return child1, child2

