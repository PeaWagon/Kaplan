"""

In charge of evolving a population for the set
number of mating events.

"""

import os
import pickle

import numpy as np

from kaplan.inputs import Inputs, InputError
from kaplan.ring import Ring, RingEmptyError
from kaplan.tournament import run_tournament
from kaplan.output import run_output
from kaplan.energy import run_energy_calc

def run_kaplan(input_dict, ring=None, save=True):
    """Run the Kaplan programme.

    Parameters
    ----------
    input_dict : dict
        Contains inputs required to run Kaplan.
    ring : str
        File name for the pickled ring object.
    save : bool
        Defaults to True, which causes the program
        to write a pickle binary file to the output
        directory containing the ring and input objects.
        False means no pickle files will be generated.
    
    Raises
    ------
    AssertionError
        The pickle file name for the ring structure
        does not exist.
        The number of slots is different than the
        original number of slots for the ring.

    """
    # read in and verify inputs
    inputs = Inputs()
    inputs.update_inputs(input_dict)

    # check that initial geometry converges
    run_energy_calc(inputs.coords)

    # biggest_bday is initial mating event number
    # it is set to the age of the oldest exsiting pmem
    biggest_bday = 0
    if ring:
        assert os.path.isfile(ring)
        assert ring.num_slots == inputs.num_slots
        for pmem in ring:
            if pmem.birthday > biggest_bday:
                biggest_bday = pmem.birthday
    else:
        # make a ring; constructor fills
        # ring with an initial population
        ring = Ring(inputs.num_slots, inputs.init_popsize)

    # run the mevs
    for mev in range(biggest_bday, biggest_bday+inputs.num_mevs):
        try:
            print(f"Mating event: {mev}")
            run_tournament(ring, mev)
            # write a temporary file every 10 mating events
            if mev%10 == 0:
                with open(os.path.join(inputs.output_dir, f"temp_{inputs.struct_input}.pickle"), "wb") as f:
                    pickle.dump(ring, f)
        except RingEmptyError:
            ring.fill(inputs.init_popsize, mev)

    # run output
    run_output(ring, save)
