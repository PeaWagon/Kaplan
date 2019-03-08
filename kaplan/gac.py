"""

In charge of evolving a population for the set
number of mating events.

In progress :D

"""

import sys

import numpy as np

from kaplan.inputs import Inputs, InputError
from kaplan.ring import Ring, RingEmptyError
from kaplan.tournament import run_tournament
from kaplan.output import run_output
from kaplan.energy import run_energy_calc

def run_kaplan(input_dict):
    """Run the Kaplan programme.

    Parameters
    ----------
    input_dict : dict
        Contains inputs required to run Kaplan.

    """
    # read in and verify inputs
    inputs = Inputs()
    inputs.update_inputs(input_dict)

    # check that initial geometry converges
    run_energy_calc(inputs.coords)

    # make a ring; constructor fills
    # ring with an initial population
    ring = Ring(inputs.num_slots, inputs.init_popsize)

    # run the mevs
    for mev in range(inputs.num_mevs):
        try:
            print(mev)
            run_tournament(ring, mev)
        except RingEmptyError:
            ring.fill(inputs.init_popsize, mev)

    # run output
    run_output(ring)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise InputError("Please include the input_dict as an argument to the program.")
    run_kaplan(sys.argv[1])
