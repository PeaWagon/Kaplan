"""
This module is responsible for generating
new population instances. The mutate function
randomly selects 0 to num_muts places to mutate
dihedral angles for each conformer geometry.
The swap takes two pmems as input and changes
the location of a random number (between 0
and num_swaps) of geometries. The generate
children function calls both mutate and
swap for two parent pmems.
NOTE:
Might want to add another mutate function
that changes the dihedral angle by +/- a set
value (to avoid massive angular changes).
"""

from random import sample, randint
from copy import deepcopy

# values for dihedral angles in degrees
MIN_VALUE = 0
MAX_VALUE = 360


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
    rng_geoms = range(len(dihedrals))
    rng_dihedrals = range(len(dihedrals[0]))
    for geom in rng_geoms:
        # choose how many mutations to do
        num_muts = randint(0, num_muts)
        # choose where to do the mutations
        mut_ind = sample(rng_dihedrals, num_muts)
        for mut in mut_ind:
            dihedrals[geom][mut] = randint(MIN_VALUE, MAX_VALUE-1)
    return dihedrals


def swap(child1, child2, num_swaps):
    """Swap geometries between two children.

    Parameters
    ----------
    child1 : list
        [[dihedrals1], [dihedrals2],
        ..., [dihedralsn]]
        Where n is the number of geometries
        in the pmem.
    child2 : list
        Same as child1, except for a different
        pmem.
    num_swaps : int
        Maximum number of swaps to do.

    Returns
    -------
    tuple of two list of lists
    Each list of list represents
    sets of dihedral angles.

    """
    # choose how many swaps to do
    num_swaps = randint(0, num_swaps)
    # choose where to do the swaps
    swap_ind = sample(range(len(child1)), num_swaps)
    # apply swaps
    for swp in swap_ind:
        child1_value = child1[swp]
        child1[swp] = child2[swp]
        child2[swp] = child1_value
    # return updated pmems
    return child1, child2
