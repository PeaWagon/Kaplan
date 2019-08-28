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
"""

from math import pi
from random import sample, randint, choice
from numpy import concatenate
from copy import deepcopy

# values for dihedral angles in radians
from kaplan.inputs import Inputs


def generate_children(parent1, parent2, num_muts, num_swaps,
                      num_cross, max_cross_points, cross_points=None):
    """Make some new pmems for the ring.

    Parameters
    ----------
    parent1 : np.array(shape=(num_geoms, num_diheds), dtype=float)
        list of lists of dihedral angles
    parent2 : np.array(shape=(num_geoms, num_diheds), dtype=float)
        list of lists of dihedral angles
    num_muts : int
        maximum number of mutations to perform
    num_swaps : int
        maximum number of swaps to perform
    num_cross : int
        maximum number of crossovers to perform.
        n-point crossover is implemented
    max_cross_points : int
        maximum number of cuts to make in
        n_point_crossover.
    cross_points : list(int)
        where to apply the n-point crossover.
        If None, crossover locations are chosen
        randomly. Each integer should be in the
        range [1, num_diheds -1] inclusive.

    Notes
    -----
    The size of parent1 and parent2 should be the same.
    num_cross and num_swaps should not exceed the length
    of parent1.

    Returns
    -------
    two new sets of dihedral angles with which
    to make pmem objects
    np.array(shape=(num_geoms, num_diheds), dtype=float)

    """
    num_conf = len(parent1)
    num_diheds = len(parent1[0])
    assert num_conf == len(parent2)
    assert num_swaps <= num_conf
    assert num_cross <= num_conf
    assert num_muts <= num_conf * num_diheds

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
        assert 0 < max_cross_points
        pmem1_indices = [i for i in range(num_conf)]
        pmem2_indices = [i for i in range(num_conf)]
        for _ in range(num_cross):
            # choose how many cuts to make
            num_cuts = randint(1, max_cross_points)
            # choose the sets of dihedrals on which to perform crossover
            conf1, conf2 = choice(pmem1_indices), choice(pmem2_indices)
            # apply crossover operator
            child1, child2 = n_point_crossover(
                child1, child2, conf1, conf2, num_cuts, cross_points
            )
            # without replacement
            pmem1_indices.remove(conf1)
            pmem2_indices.remove(conf2)

    if num_muts1:
        # comb is all the possible locations that can be mutated in the child pmem
        comb = [(i, j) for i in range(num_conf) for j in range(num_diheds)]
        # select a random location num_muts times
        mutate_locations = sample(comb, num_muts1)
        for loc in mutate_locations:
            # apply mutate operator to child1
            child1 = mutate(child1, loc)

    if num_muts2:
        comb = [(i, j) for i in range(num_conf) for j in range(num_diheds)]
        mutate_locations = sample(comb, num_muts2)
        for loc in mutate_locations:
            child2 = mutate(child2, loc)

    return child1, child2


def single_parent_mutation(parent, num_muts):
    """Mutate a parent.

    Parameters
    ----------
    parent : np.array(shape=(num_geoms, num_diheds), dtype=float)
        list of lists of dihedral angles
    num_muts : int
        maximum number of mutations to perform

    Returns
    -------
    parent with x point mutations, where x <= num_muts

    """
    num_conf = len(parent)
    num_diheds = len(parent[0])
    child = deepcopy(parent)
    # choose how many mutations to do
    num_muts = randint(0, num_muts)
    if num_muts:
        # comb is all the possible locations that can be mutated in the child pmem
        comb = [(i, j) for i in range(num_conf) for j in range(num_diheds)]
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
    inputs = Inputs()
    geom, dihed = loc
    dihedrals[geom][dihed] = choice(inputs.avail_diheds)
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


def n_point_crossover(child1, child2, conf1, conf2, n, cut_indices=None):
    """Perform single point crossover n times.

    Parameters
    ----------
    child1 : np.array(shape=(num_geoms, num_diheds), dtype=float)
        Set of dihedral angles for each conformer in child1.
    child2 : np.array(shape=(num_geoms, num_diheds), dtype=float)
        Set of dihedral angles for each conformer in child2.
    conf1 : int
        Location of conformer 1 in child1 for crossover.
    conf2 : int
        Location of conformer 2 in child2 for crossover.
    n : int
        How many places to cut the set of dihedral angles.
    cut_indices : list(int)
        Defaults to None, which means the cut indices are
        generated by randomly sampling (without replacement)
        from the half interval [1, num_diheds). If cut
        indices are given, they should be integer values
        in the same half interval [1, num_diheds).

    Algorithm
    ---------
    Randomly select conformer ca from child1.
    Randomly select conformer cb from child2.
    Repeat n times:
        Randomly select index d for crossover, such that
        the separated subarrays cannot be empty (i.e. index
        cannot be 0 or len(array)).
        ca is now ca[:d] + cb[d:]
        cb is now cb[:d] + ca[d:]
    Return 2 new children.

    Example
    -------
    Here is an example of single point crossover (n == 1):
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

    Notes
    -----
    If the num_diheds is 1, then crossover cannot be completed.

    Returns
    -------
    2 x np.array(shape=(num_geoms, num_diheds), dtype=float)

    """
    # n is how many cuts to make
    num_diheds = len(child1[0])
    # uniform crossover would be the case where this is equal
    assert 0 < n <= num_diheds - 1
    if num_diheds == 1:
        print(f"Warning. {n}-point crossover could not be completed.")
        print("If the molecule only has one dihedral angle, try swap instead.")
        return child1, child2
    if cut_indices is None:
        cut_indices = sample(range(1, num_diheds), n)
    else:
        # cuts at the same location are allowed, but they will simply
        # cause a null operation (no change)
        cut_indices = sample(cut_indices, n)
        # make sure cuts are within allowed range
        for cut in cut_indices:
            assert isinstance(cut, int)
            assert cut in range(1, num_diheds)

    for cut in cut_indices:
        # create two new sets of dihedrals to put back into the child pmems
        new1 = deepcopy(child1[conf1])
        new2 = deepcopy(child2[conf2])
        # with numpy, adding two arrays adds element-wise
        # need to use concatenate to achieve list-like addition
        child1[conf1] = concatenate([new1[:cut], new2[cut:]])
        child2[conf2] = concatenate([new2[:cut], new1[cut:]])

    return child1, child2
