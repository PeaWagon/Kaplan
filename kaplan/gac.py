"""

In charge of evolving a population for the set
number of mating events.

"""

import sys

import numpy as np

from kaplan.inputs import Inputs, InputError
from kaplan.ring import Ring, RingEmptyError
from kaplan.tournament import run_tournament
from kaplan.output import run_output
from kaplan.energy import run_energy_calc

def run_kaplan(input_dict, save=True, output_loc="pwd"):
    """Run the Kaplan programme.

    Parameters
    ----------
    input_dict : dict
        Contains inputs required to run Kaplan.
    save : bool
        Defaults to True, which causes the program
        to write a pickle binary file to the output
        directory containing the ring and input objects.
        False means no pickle files will be generated.
    output_loc : str
        Where to save the output. Defaults to present/
        current working directory, under kaplan_output
        followed by a job directory. If home is specified,
        the output for the first job is written to:
        /home/username/kaplan_output/job_0).

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
            print(f"Mating event: {mev}")
            run_tournament(ring, mev)
        except RingEmptyError:
            ring.fill(inputs.init_popsize, mev)

    # run output
    run_output(ring, save, output_loc)
