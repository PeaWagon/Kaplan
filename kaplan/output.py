
OUTPUT_FORMAT = 'xyz'

def run_output(ring):
    """Run the output module.

    Parameters
    ----------
    ring : object
       the final ring data structure after evolution.

    """
    # find average fitness
    # find best pmem and its ring index
    total_fit = 0
    best_pmem = 0
    best_fit = 0
    for pmem in ring.pmems:
        if pmem is not None:
            fitg = pmem.fitness
            total += fitg
            if fitg > best_fit:
                best_pmem = pmem.ring_loc
                best_fit = fitg
    average_fit = total_fit / ring.num_filled

    # generate the output file for the best pmem


