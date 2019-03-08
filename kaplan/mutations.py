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

from math import pi
from random import sample, randint, choice
from numpy.random import uniform
from numpy import concatenate
from copy import deepcopy

# values for dihedral angles in radians
MIN_VALUE = 0
MAX_VALUE = 2*pi


def generate_children(parent1, parent2, num_muts, num_swaps, cross=False):
    """Make some new pmems for the ring.

    Parameters
    ----------
    parent1 : np.array(shape=(num_geoms, num_dihed), dtype=float)
        list of lists of dihedral angles
    parent2 : np.array(shape=(num_geoms, num_dihed), dtype=float)
        list of lists of dihedral angles
    num_muts : int
        maximum number of mutations to perform
    num_swaps : int
        maximum number of swaps to perform
    cross : bool
        True means do single point crossover. False
        means don't do single point crossover.

    Returns
    -------
    two new sets of dihedral angles with which
    to make pmem objects
    np.array(shape=(num_geoms, num_dihed), dtype=float)

    """
    child1 = deepcopy(parent1)
    child2 = deepcopy(parent2)
    child1, child2 = swap(child1, child2, num_swaps)
    if cross:
        child1, child2 = crossover(child1, child2)
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
            dihedrals[geom][mut] = uniform(MIN_VALUE, MAX_VALUE)
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


def crossover(child1, child2):
    """Perform single point crossover between two children.

    Parameters
    ----------
    child1 : np.array(shape=(num_geoms, num_dihed), dtype=float)
        Set of dihedral angles for each conformer in child1.
    child2 : np.array(shape=(num_geoms, num_dihed), dtype=float)
        Set of dihedral angles for each conformer in child2.
    The size of child1 and child2 should be the same.

    Algorithm
    ---------
    Randomly select conformer ca from child1.
    Randomly select conformer cb from child2.
    Randomly select index d for crossover, such that
    the separated subarrays cannot be empty (i.e. index
    cannot be 0 or len(array)).
    ca is now ca[:d] + cb[d:]
    cb is now cb[:d] + ca[d:]
    Return 2 new children.

    Example
    -------
    child1 = np.array([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]])
    child2 = np.array([[13, 14, 15, 16], [17, 18, 19, 20], [21, 22, 23, 24]])
    cut index: 3
    chosen1: 2
    chosen2: 1
    new1: [ 9 10 11 12] 
    new2: [17 18 19 20]
    child1: [[ 1  2  3  4]
             [ 5  6  7  8]
             [ 9 10 11 20]]
    child2: [[13 14 15 16]
             [17 18 19 12]
             [21 22 23 24]]
    
    Returns
    -------
    2 x np.array(shape=(num_geoms, num_dihed), dtype=float)

    """
    cut_index = randint(1, (len(child1[0]) - 1))
    max_select = len(child1) - 1
    chosen1 = randint(0, max_select)
    chosen2 = randint(0, max_select)
    new1 = deepcopy(child1[chosen1])
    new2 = deepcopy(child2[chosen2])
    # with numpy, adding two arrays adds element-wise
    # need to use concatenate to achieve list-like addition
    child1[chosen1] = concatenate([new1[:cut_index], new2[cut_index:]])
    child2[chosen2] = concatenate([new2[:cut_index], new1[cut_index:]])
    return child1, child2
