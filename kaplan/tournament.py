
import numpy as np
from kaplan.ring import RingEmptyError
from kaplan.mutations import generate_children

def run_tournament(t_size, num_muts, num_swaps, ring,
                   current_mev):
    # check ring has enough pmems for a tournament
    if t_size > ring.num_filled:
        raise RingEmptyError("Not enough pmems to run a tournament.")

    # choose random slots for a tournament
    selected_pmems = select_pmems(t_size, ring)

    # select parents by fitness
    parents = select_parents(selected_pmems, ring)
    print(parents)
    print(parents[0])
    print(ring.pmems)
    print('h')
    print(ring.pmems[parents[0]].dihedrals)
    print('ej')
#    parent1 = ring.pmems[next(parents)].dihedrals
#    parent2 = ring.pmems[next(parents)].dihedrals
    parent1 = ring.pmems[parents[0]].dihedrals
    parent2 = ring.pmems[parents[1]].dihedrals

    # generate children
    children = generate_children(parent1, parent2, num_muts,
                                 num_swaps)

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

    Note
    ----
    Maybe move this to the ring as a method. Or turn
    into a generator (calling next on the ring returns
    a random pmem).

    """
    selection = []
    for pmem in range(number):
        # choose random slot
        choice = np.random.randint(0, ring.num_slots)
        # if the slot is empty, keep picking another random one
        while ring.pmems[choice] is None:
            choice = np.random.randint(0, ring.num_slots)
        selection.append(choice)
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
    tuple : int,int
        Indices of two best pmems to be used as parents.

    """
    fit_vals = np.array([ring.pmems[i].fitness for i in selected_pmems])
    print(fit_vals)
    # from here:
    # https://stackoverflow.com/questions/6910641/how-do-i-get-indices-of-n-maximum-values-in-a-numpy-array
    # use numpy to get the two best fitness value indices
    # from the list and link it to ring index
    # note: this makes generator object
    parents_gen = (selected_pmems[parent] for parent in np.argpartition(fit_vals, -2)[-2:])
    parents = [next(parents_gen)]
    parents.append(next(parents_gen))
    print(parents)
    return parents

