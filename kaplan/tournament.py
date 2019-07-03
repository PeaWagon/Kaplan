"""The tournament module decides which pmems to
pick from the ring in order to apply updates
to the population."""

from random import choice

from kaplan.inputs import Inputs, InputError
from kaplan.ring import RingEmptyError
from kaplan.mutations import generate_children, single_parent_mutation


def run_tournament(ring, current_mev, test=False):
    """Run the tournament (i.e. a mating event).

    Parameters
    ----------
    ring : object
        Ring object.
    current_mev : int
        The current mating event number. Used
        to give pmems birthdays.
    test : bool
        Defaults to False. If set to True, extra information
        is printed to the terminal detailing the tournament
        participants and results.

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
    selected_pmems, parent1_loc = select_pmems(inputs.mating_rad, ring)

    if test:
        print("Pmems in tournament:", selected_pmems)
        print("Parent 1:", parent1_loc)

    # if normalise is selected, apply it to the tournament selection
    # the fitness of parent1 does not need to be updated at this time
    # it should only be updated if the parent1's slot is available
    # as a slot for child1 or child2 (right now, the minimum number
    # of slots for a tournament is 5 - the user can force this to
    # change, but it would probably break the program)
    if inputs.normalise:
        for slot in selected_pmems:
            if ring[slot] is not None:
                ring.set_fitness(ring[slot])
    p1_dihed = ring[parent1_loc].dihedrals

    # select parents and worst slots by fitness
    parent2_loc, worst = select_parents(selected_pmems, ring)

    if test:
        print("Worst pmem locations:", worst)
        print("Parent 2:", parent2_loc)

    try:
        p2_dihed = ring[parent2_loc].dihedrals
    # there is only one pmem in the mating radius
    except AttributeError:
        child1 = single_parent_mutation(p1_dihed, inputs.num_muts)
        child2 = None
    else:
        # generate children
        child1, child2 = generate_children(
            p1_dihed, p2_dihed,
            inputs.num_muts, inputs.num_swaps,
            inputs.num_cross
        )
    # put children in ring
    ring.update(child1, worst[0], current_mev)
    if child2 is not None:
        ring.update(child2, worst[1], current_mev)


def select_pmems(mating_rad, ring):
    """Randomly selected pmems.

    Parameters
    ----------
    mating_rad : int
        How many slots (to the left and right) to
        consider when selecting pmems from the ring.
    ring : object
        The ring from which to pick the pmems.

    Returns
    -------
    tuple : (list(int), int)
        List of ring indices with which to run the tournament and
        the index of parent1.

    """
    # get a list of indices representing filled slots
    occupied = [i for i in range(ring.num_slots) if ring[i] is not None]
    # from the occupied slots, choose a random pmem as parent1
    parent1 = choice(occupied)
    selection = ring.mating_radius(parent1, mating_rad)
    # don't include parent1 in the selection, so its fitness is not
    # compared to the other pmems
    # this choice gives low-fitness pmems a chance to mate
    selection.remove(parent1)
    return selection, parent1


def sortby(ring, sortkey, include_slots="all"):
    """Sort a selection of ring slots by pmem attribute.

    Parameters
    ----------
    ring : kaplan.ring.Ring object
        Contains pmems that should be sorted.
    sortkey : str
        How to sort the pmems in the ring. The sortkey
        should be a valid pmem attribute. Examples
        include fitness and birthday.
    include_slots : list(int)
        Each integer is a ring index that is included
        in the sorting. Defaults to "all", which means
        all slots are included.

    Raises
    ------
    AssertionError
        The number of slots to include exceeds the
        ring's capacity (in terms of number of slots).
        Also raised if there is a slot number that
        is not within the ring's available slots.

    Returns
    -------
    index_val_pairs : list(tuple(int, float or None))
        A list containing sorted slot indices, with
        the value given as a float or None (for empty
        or pmems with uninitialised values).
        The slots start lowest values (or None)
        to highest values.

    """
    if include_slots == "all":
        include_slots = range(ring.num_slots)
    else:
        # remove duplicated slot numbers
        include_slots = set(include_slots)
        # make sure slots are valid
        assert len(include_slots) <= ring.num_slots
        assert all(slot in range(ring.num_slots) for slot in include_slots)

    index_val_pairs = []
    empty_slots = []
    last_index = -1
    # values are tuple(ring index, fitness)
    for s in include_slots:
        # attribute error should occur here if invalid
        # sortkey is called
        if ring[s] is None or getattr(ring[s], sortkey) is None:
            empty_slots.append((s, None))
        else:
            last_index += 1
            index_val_pairs.append((s, getattr(ring[s], sortkey)))
    quicksort(index_val_pairs, 0, last_index)
    index_val_pairs = empty_slots + index_val_pairs
    return index_val_pairs


def select_parents(selected_pmems, ring):
    """Sort pmems by fitness values, choose best 2.

    Parameters
    ----------
    selected_pmems : list
        The indices of the ring that are in the tournament.
        Not including parent 1.

    ring : object
        Instance of Ring class.

    Notes
    -----
    If a pmem does not have its fitness set (i.e. fitness
    value of None), it will be treated as an empty slot
    and placed at the worst end of the list during fitness
    sorting.

    Returns
    -------
    tuple : int, tuple(int, int)
        Parent2 index and indices of two worst pmems
        to potentially be replaced.

    """
    index_fit_pairs = sortby(ring, "fitness", selected_pmems)
    worst_two = index_fit_pairs[0][0], index_fit_pairs[1][0]
    parent2 = index_fit_pairs[-1][0]
    return parent2, worst_two


def quicksort(pairs, low, high):
    """Sort a list of index-value pairs.

    Parameters
    ----------
    pairs : list(tuple(int, float))
        each pair has the following format: (location, value)
        where location is slot number in the ring and value is fitness
    low : int
        Lowest index to consider from pairs to sort.
    high : int
        Highest index to consider from pairs to sort.
        Should NOT be equal to the length of pairs.

    Notes
    -----
    This is the quicksort algorithm that starts by setting the pivot
    value to the last item in the list. It is an in-place sorting,
    meaning that it does not return a value, but instead changes
    the order of pairs. It is a recursive algorithm.

    """
    if low < high:
        # after partition, pi is in the correct place
        pi = partition(pairs, low, high)
        quicksort(pairs, low, pi - 1)
        quicksort(pairs, pi + 1, high)


def partition(pairs, low, high):
    """Partition the pairs from low to high relative to the pivot.

    Parameters
    ----------
    Same as quicksort.

    Notes
    -----
    Any item in pairs prior to the pivot value that is lower
    in value to the pivot is placed at the start of pairs.
    How many items that are less in value than the pivot are
    counted. Then the pivot is placed after the number of items
    found smaller in value to it.

    Returns
    -------
    i : int
        The index of the correctly placed item (after partition).

    """
    pivot = pairs[high]
    i = (low - 1)
    for j in range(low, high):
        if pairs[j][1] <= pivot[1]:
            i += 1
            swap(pairs, i, j)
    swap(pairs, i + 1, high)
    return i + 1


def swap(pairs, i, j):
    """Given pairs, swap locations i and j in the list."""
    valuei = pairs[i]
    pairs[i] = pairs[j]
    pairs[j] = valuei
