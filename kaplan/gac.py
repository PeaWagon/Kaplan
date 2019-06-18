"""

In charge of evolving a population for the set
number of mating events.

"""

import os
import pickle

import numpy as np

from kaplan.inputs import Inputs, InputError, read_input
from kaplan.ring import Ring, RingEmptyError
from kaplan.tournament import run_tournament
from kaplan.output import run_output
from kaplan.energy import run_energy_calc

def run_kaplan(inputs, ring=None, save=True, new_dir=True):
    """Run the Kaplan programme.

    Parameters
    ----------
    inputs : dict or str
        Contains inputs required to run Kaplan.
        If dict, key/value pairs are used to generate
        the inputs object. If str, then this should
        be the path and file name of the pickled inputs
        object to read as input (i.e. inputs.pickle).
    ring : str
        Path and file name for the pickled ring object
        (i.e. ring.pickle).
    save : bool
        Defaults to True, which causes the program
        to write a pickle binary file to the output
        directory containing the ring and input objects.
        False means no pickle files will be generated.
    new_dir : bool
        Defaults to True, which means that a new job
        directory will be generated in the case that
        an old inputs and/or ring object is used. If
        False, then the original directory will be
        kept (if inputs pickle is given as input).

    
    Raises
    ------
    AssertionError
        The pickle file name for the ring structure
        does not exist.
        The number of slots is different than the
        original number of slots for the ring.

    """
    # read in and verify inputs
    if isinstance(inputs, str):
        assert os.path.isfile(inputs)
        inputs = read_input(inputs, new_dir)
    else:
        inputs = Inputs()
        inputs.update_inputs(input_dict)

    # check that initial geometry converges
    result = run_energy_calc(inputs.coords)
    if result == -1:
        print("Warning: initial geometry did not converge.")

    # biggest_bday is initial mating event number
    # it is set to the age of the oldest exsiting pmem
    biggest_bday = 0
    if ring:
        assert os.path.isfile(ring)
        with open(ring, "rb") as f:
            ring = pickle.load(f)
        assert ring.num_slots == inputs.num_slots
        for pmem in ring:
            if pmem.birthday > biggest_bday:
                biggest_bday = pmem.birthday
    else:
        # make a ring; constructor fills
        # ring with an initial population
        ring = Ring(inputs.num_slots, inputs.init_popsize)

    temp_name = inputs.struct_input
    # strip the path and extension from the struct_input
    if inputs.struct_type in ("xyz", "com"):
        temp_name = os.path.basename(os.path.splitext(temp_name)[0])

    # run the mevs
    for mev in range(biggest_bday, biggest_bday+inputs.num_mevs):
        try:
            print(f"Mating event: {mev}")
            run_tournament(ring, mev)
            # write a temporary file every 10 mating events
            if mev%10 == 0:
                with open(os.path.join(inputs.output_dir, f"temp_{temp_name}.pickle"), "wb") as f:
                    pickle.dump(ring, f)
        except RingEmptyError:
            ring.fill(inputs.init_popsize, mev)

    # run output
    run_output(ring, save)
