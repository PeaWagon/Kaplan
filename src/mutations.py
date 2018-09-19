

import numpy as np

# values for dihedral angles in degrees
MIN_VALUE = 0
MAX_VALUE = 360

"""
Might want to add another mutate function
that changes the dihedral angle by +/- a set
value (to avoid massive angular changes).

"""

def mutate(input_list, min_muts, max_muts):
    # determine explicit number of mutations
    num_muts = np.random.randint(min_muts, max_muts+1, size=1)
    # if num_muts is larger than input_list length,
    # do mutation with replacement, otherwise do
    # mutation without replacement
    if len(input_list) > num_muts:
        # determine indices of mutation
        mut_indices = []
        for m in mut_indices:
            # generate new value (can be equal to prev value)
            new = []
            # replace old value
            input_list[m] = new
        # return new mutated list
        return input_list
    # generate entirely new input_list; return
    else:
        return []

def crossover(list1, list2, num_cuts):
    pass

