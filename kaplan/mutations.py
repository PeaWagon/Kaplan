

import numpy as np
from random import sample

# values for dihedral angles in degrees
MIN_VALUE = 0
MAX_VALUE = 360

"""
Might want to add another mutate function
that changes the dihedral angle by +/- a set
value (to avoid massive angular changes).

"""

def mutate(input_list, min_muts, max_muts, test=False):
    """Mutate an input list of dihedral angles.

    Parameters
    ----------
    input_list : list
        A list of dihedral angles (integers).
    min_muts : int
        The minimum number of mutations to perform.
    max_muts : int
        The maximum number of mutations to perform.
    test : bool, list
        If False, no test is being performed.
        Otherwise, test is a list of lists, where
        the first value is a list specifying the
        locations to mutate and the second value
        is a list of new values to use in the
        returned list.

    Returns
    -------
    A list of dihedral angles after mutation.

    Notes
    -----
    If the number of mutations is greater than the
    number of elements in the input_list, then the
    whole list will be regenerated. The number of
    mutations is chosen randomly, starting with
    min_muts up until max_muts.

    The new values for the dihedral angles will be
    any integer value between the minimum and the
    maximum value (as specified by the global
    variables MIN_VALUE and MAX_VALUE).

    Raises
    ------
    AssertionError
        The max_muts value is smaller than the
        min_muts value.

    """
    assert max_muts >= min_muts
    if not test:
        # determine explicit number of mutations
        num_muts = np.random.randint(min_muts, max_muts+1)
    else:
        # make sure test is not doing something stupid
        assert len(test[0]) == len(test[1])
        num_muts = len(test[0])
    # if num_muts is larger than input_list length,
    # do mutation with replacement, otherwise do
    # mutation without replacement
    if len(input_list) > num_muts or test:
        # determine indices of mutation
        if not test:
            mut_indices = sample(range(len(input_list)), num_muts)
        else:
            mut_indices = test[0]
        for i, m in enumerate(mut_indices):
            if not test:
                # new value can be equal to prev value
                new = np.random.randint(MIN_VALUE, MAX_VALUE)
            else:
                new = test[1][i]
            # replace old value
            input_list[m] = new
        # return new mutated list
        return input_list
    # generate entirely new input_list; return
    else:
        return sample(range(MIN_VALUE, MAX_VALUE), len(input_list))

def crossover(list1, list2, num_cuts):
    min_len = min(len(list1), len(list2))
    assert num_cuts <= min_len - 1















