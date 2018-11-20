
import numpy as np

def tournament():
    pass

def set_tourney(t_size, ring):
    """Determine which pmems are seleccted from the ring."""
    # check ring has enough pmems for a tournament
    # if ring.num_filled >= t_size:
    try:
        # Generate a uniform random sample from num_filled of
        # size t_size without replacement, aka don't pick
        # same pmem twice
        selected_pmems = np.random.choice(ring.num_filled,
                                          t_size, replace=False)
        return selected_pmems
    # Cannot take a larger sample than population when
    # 'replace=False'
    # aka size of set is larger than available pmems
    except ValueError:
        # generate new pmems randomly or by using available pmems
        pass
    


def sort_fitg(selected_pmems, ring):
    
