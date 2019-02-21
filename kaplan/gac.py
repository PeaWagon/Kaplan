"""

In charge of evolving a population for the set
number of mating events.

In progress :D

"""

import sys

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

    # keep track of initial population size
    # in case ring needs to be refilled
    init_popsize = inputs.num_filled

    # check that initial geometry converges
    assert run_energy_calc(inputs.parser.coords)

    # make a ring
    ring = Ring()

    # fill ring with an initial population
    ring.fill(init_popsize, 0)

    # run the mevs
    for mev in range(inputs.num_mevs):
        try:
            print(mev)
            run_tournament(ring, mev)
        except RingEmptyError:
            ring.fill(init_popsize, mev)

    # run output
    run_output(ring)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise InputError("Please include the input_dict as an argument to the program.")
    run_kaplan(sys.argv[1])
