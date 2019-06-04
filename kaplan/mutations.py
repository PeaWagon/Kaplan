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
from kaplan.geometry import MAX_VALUE, MIN_VALUE



def generate_children(parent1, parent2, num_muts, num_swaps, num_cross):
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
    num_cross : int
        maximum number of crossovers to perform.
        Single-point crossover is implemented
    
    Notes
    -----
    The size of parent1 and parent2 should be the same.
    num_cross and num_swaps should not exceed the length
    of parent1.
    
    Returns
    -------
    two new sets of dihedral angles with which
    to make pmem objects
    np.array(shape=(num_geoms, num_dihed), dtype=float)

    """
    num_conf = len(parent1)
    num_dihed = len(parent1[0])
    assert num_conf == len(parent2)
    assert num_swaps <= num_conf
    assert num_cross <= num_conf
    assert num_muts <= num_conf*num_dihed
    
    child1 = deepcopy(parent1)
    child2 = deepcopy(parent2)
    # choose how many swaps to do
    num_swaps = randint(0, num_swaps)
    # choose how many crossovers to do
    num_cross = randint(0, num_cross)
    # choose how many mutations to do
    num_muts1 = randint(0, num_muts)
    num_muts2 = randint(0, num_muts)

    if num_swaps:
        pmem1_indices = [i for i in range(num_conf)]
        pmem2_indices = [i for i in range(num_conf)]
        for _ in range(num_swaps):
            # choose where to perform swap
            conf1, conf2 = choice(pmem1_indices), choice(pmem2_indices)
            # apply swap operator
            child1, child2 = swap(child1, child2, conf1, conf2)
            # without replacement
            pmem1_indices.remove(conf1)
            pmem2_indices.remove(conf2)
        
    if num_cross:
        pmem1_indices = [i for i in range(num_conf)]
        pmem2_indices = [i for i in range(num_conf)]
        for _ in range(num_cross):
            # choose where to perform crossover
            conf1, conf2 = choice(pmem1_indices), choice(pmem2_indices)
            # apply crossover operator
            child1, child2 = crossover(child1, child2, conf1, conf2)
            # without replacement
            pmem1_indices.remove(conf1)
            pmem2_indices.remove(conf2)

    if num_muts1:
        # comb is all the possible locations that can be mutated in the child pmem
        comb = [(i,j) for i in range(num_conf) for j in range(num_dihed)]
        # select a random location num_muts times
        mutate_locations = sample(comb, num_muts1)
        for loc in mutate_locations:
            # apply mutate operator to child1
            child1 = mutate(child1, loc)
    
    if num_muts2:
        comb = [(i,j) for i in range(num_conf) for j in range(num_dihed)]
        mutate_locations = sample(comb, num_muts2)
        for loc in mutate_locations:
            child2 = mutate(child2, loc)
    
    return child1, child2


def single_parent_mutation(parent, num_muts):
    """Mutate a parent.

    Parameters
    ----------
    parent : np.array(shape=(num_geoms, num_dihed), dtype=float)
        list of lists of dihedral angles
    num_muts : int
        maximum number of mutations to perform

    Returns
    -------
    parent with x point mutations, where x <= num_muts

    """
    num_conf = len(parent)
    num_dihed = len(parent[0])
    child = deepcopy(parent)
    # choose how many mutations to do
    num_muts = randint(0, num_muts)
    if num_muts:
        # comb is all the possible locations that can be mutated in the child pmem
        comb = [(i,j) for i in range(num_conf) for j in range(num_dihed)]
        # select a random location num_muts times
        mutate_locations = sample(comb, num_muts)
        for loc in mutate_locations:
            # apply mutate operator to child1
            child = mutate(child, loc)
    return child


def mutate(dihedrals, loc):
    """Mutate a list of dihedral angles.

    Parameters
    ----------
    dihedrals : list(list (int))
        List of lists of dihedral angles to mutate.
    loc : tuple(int, int)
        Where in the dihedrals to mutate.

    Returns
    -------
    dihedrals : list(list (int))
        Dihedral angles after mutations.

    """
    geom, dihed = loc
    dihedrals[geom][dihed] = uniform(MIN_VALUE, MAX_VALUE)
    return dihedrals


def swap(child1, child2, conf1, conf2):
    """Swap a conformer between two children.

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
    conf1 : int
        Location of conformer 1 in child1 to be
        swapped.
    conf2 : int
        Location of conformer 2 in child2 to be
        swapped.

    Returns
    -------
    tuple of two list of lists
    Each list of list represents
    sets of dihedral angles.

    """  
    child1_value = deepcopy(child1[conf1])
    child1[conf1] = child2[conf2]
    child2[conf2] = child1_value
    # return updated pmems
    return child1, child2


def crossover(child1, child2, conf1, conf2):
    """Perform single point crossover for conformers from two children.

    Parameters
    ----------
    child1 : np.array(shape=(num_geoms, num_dihed), dtype=float)
        Set of dihedral angles for each conformer in child1.
    child2 : np.array(shape=(num_geoms, num_dihed), dtype=float)
        Set of dihedral angles for each conformer in child2.
    conf1 : int
        Location of conformer 1 in child1 for crossover.
    conf2 : int
        Location of conformer 2 in child2 for crossover.

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
    new1 = deepcopy(child1[conf1])
    new2 = deepcopy(child2[conf2])
    # with numpy, adding two arrays adds element-wise
    # need to use concatenate to achieve list-like addition
    child1[conf1] = concatenate([new1[:cut_index], new2[cut_index:]])
    child2[conf2] = concatenate([new2[:cut_index], new1[cut_index:]])
    return child1, child2
