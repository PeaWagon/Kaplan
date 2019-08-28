"""

In charge of evolving a population for the set
number of mating events.

"""

import os
import pickle

from shutil import copyfile
from random import random
from copy import deepcopy

from numpy import allclose

from kaplan.inputs import Inputs, InputError, read_input
from kaplan.ring import Ring, RingEmptyError
from kaplan.tournament import run_tournament
from kaplan.extinction import asteroid, agathic, deluge, plague
from kaplan.output import run_output
from kaplan.energy import run_energy_calc
from kaplan.tools import make_2d, plot_2d, plot_3d
from kaplan.optimise import optimise_coords
from kaplan.geometry import write_coords
from kaplan.fitness import update_all_fitness


def run_kaplan(job_inputs, ring=None, save=True, new_dir=True,
               save_every=25, recalc_fit=True):
    """Run the Kaplan programme.

    Parameters
    ----------
    job_inputs : dict or str or kaplan.inputs.Inputs object
        Contains inputs required to run Kaplan.
        If dict, key/value pairs are used to generate
        the inputs object. If str, then this should
        be the path and file name of the pickled inputs
        object to read as input (i.e. inputs.pickle).
        If inputs object, the object's _check_inputs
        method is called prior to execution.
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
    save_every : int
        After x number of mating events, the stats_file.txt
        will be updated and a temporary ring object will
        be written to output_dir (if save is True).
    recalc_fit : bool
        This input parameter only applies if a ring object
        is input (i.e. ring is not None). If True (which
        is the default), each pmem in the ring has
        its fitness recalculated. This also will update
        the conformer energies for each pmem. The RMSD
        should not change, regardless of input. The cases
        in which this update is necessary are those where
        the input parameters change related to fitness.
        If False, the pmem fitness values will not be
        re-evaluated.

    Notes
    -----
    recalc_fit should be set to True if re-running a ring
    object after changing one or more of the following
    inputs:
    * prog/method/basis
    * coef_energy/coef_rmsd
    * fit_form
    * normalise
    * anything to do with pmem optimisation

    Raises
    ------
    AssertionError
        The pickle file name for the ring structure
        does not exist.
        The number of slots is different than the
        original number of slots for the ring.

    """
    # make sure no divide by zero errors
    assert save_every > 0

    # read in and verify inputs
    if isinstance(job_inputs, str):
        assert os.path.isfile(job_inputs)
        inputs = read_input(job_inputs, new_dir)
    elif isinstance(job_inputs, Inputs):
        inputs = job_inputs
        if new_dir:
            inputs = read_input(inputs, new_dir)
        inputs._check_input()
    else:
        inputs = Inputs()
        inputs.update_inputs(job_inputs)

    # if requested (default), optimise the initial geometry before
    # setting the initial coordinates
    # check that initial geometry converges
    # note error is not handled here, so any errors
    # will stop program execution
    if inputs.opt_init_geom:
        outfile = os.path.join(inputs.output_dir, "input_coords.xyz")
        # if it exists, make a copy of the original input_coords
        if os.path.isfile(outfile):
            copyfile(outfile, os.path.join(inputs.output_dir, "old_input_coords.xyz"))
        print("Optimising initial geometry.")
        _, coords = optimise_coords(inputs.coords, major=True)
        # update the input coordinates, and write over the
        # input coordinates xyz file
        inputs.coords = coords
        comments = "input coords after optimisation"
        write_coords(coords, inputs.atomic_nums, outfile, comments)
    else:
        run_energy_calc(inputs.coords)

    # generate plots for the input coordinates
    # in the output_dir
    make_2d()
    plot_2d()
    plot_3d()

    # biggest_bday is initial mating event number
    # it is set to the age of the oldest exsiting pmem, plus one
    biggest_bday = 0
    if ring:
        assert os.path.isfile(ring)
        with open(ring, "rb") as f:
            ring = pickle.load(f)
        assert ring.num_slots == inputs.num_slots
        for pmem in ring:
            if pmem.birthday > biggest_bday:
                biggest_bday = pmem.birthday
            if recalc_fit:
                pmem.setup(major=False)
        biggest_bday += 1
    else:
        # make a ring; constructor fills
        # ring with an initial population
        ring = Ring(inputs.num_slots, inputs.init_popsize)

    def mating_event(mev):
        """Run a mating event.

        Parameters
        ---------
        mev : int
            The mating event number to run.

        Returns
        -------
        needs_fitness_update : bool
            True if fitness needs updating and was not
            updated during the mating event. If False,
            fitness was already updated and does not
            need to be rerun.

        """
        # keep track of whether fitness values
        # need to be set
        needs_fitness_update = True
        try:
            print(f"Mating event: {mev}")
            run_tournament(ring, mev)

            # extinction events
            if inputs.asteroid:
                percent_chance = random()
                if inputs.asteroid >= percent_chance:
                    print("Applying asteroid")
                    asteroid(ring, ring.occupied)
            if inputs.plague:
                percent_chance = random()
                if inputs.plague >= percent_chance:
                    print("Applying plague")
                    # first update fitness values
                    update_all_fitness(ring)
                    needs_fitness_update = False
                    plague(ring)
            if inputs.agathic:
                percent_chance = random()
                if inputs.agathic >= percent_chance:
                    print("Applying agathic")
                    agathic(ring)
            if inputs.deluge:
                percent_chance = random()
                if inputs.deluge >= percent_chance:
                    print("Applying deluge")
                    # first update fitness values
                    update_all_fitness(ring)
                    needs_fitness_update = False
                    deluge(ring, ring.occupied)

            # keep data every 25 mevs
            if mev % save_every == 0:
                # update fitness values as required
                # normalised fitness will change if extinction
                # events deluge or plague occurred, since the population
                # size will change and thus change the mean/stdev values
                if needs_fitness_update or inputs.normalise:
                    update_all_fitness(ring)
                    needs_fitness_update = False
                run_output(ring, mev, save)

        except RingEmptyError:
            ring.fill(inputs.init_popsize, mev)

        return needs_fitness_update

    # determine if we are running until convergence
    # or a set number of mating events
    if not inputs.stop_at_conv:
        for mev in range(biggest_bday, biggest_bday + inputs.num_mevs):
            mating_event(mev)

    else:
        mev = biggest_bday
        # properties of the best pmem that was last checked
        last_best_energies = None
        last_best_rmsds = None
        while True:
            needs_fitness_update = mating_event(mev)
            # see if the energies and RMSD values of the best
            # pmem are the same
            # the RMSD and energies are sorted to ensure that
            # the geometry order does not influence the convergence
            # i.e. if two pmems are identical except for order
            # of geometries, energies and RMSD values will be
            # the same but in a different order
            if mev % inputs.stop_at_conv == 0 and mev != 0:
                if needs_fitness_update:
                    update_all_fitness(ring)
                # find best pmem and its ring index
                # this also updates all the fitness values
                # for all pmems in the ring
                best_pmem, _ = ring.best_pmem

                # if any are None, we assume that the algorithm
                # should still be running (since the geometries
                # are not all viable)
                if any(e is None for e in ring[best_pmem].energies):
                    continue
                if inputs.num_geoms != 1:
                    if any(r[2] is None for r in ring[best_pmem].rmsds):
                        continue

                # get the energies and RMSD values from the best pmem
                # as a deepcopy so we don't modify the values
                # or order in the ring
                best_energies = deepcopy(ring[best_pmem].energies)
                best_rmsds = deepcopy([rmsd[2] for rmsd in ring[best_pmem].rmsds])

                best_energies = sorted(best_energies)
                best_rmsds = sorted(best_rmsds)
                try:
                    if allclose(last_best_energies, best_energies) and\
                       allclose(last_best_rmsds, best_rmsds):
                        break
                # if last_best_x is None
                # TypeError: ufunc 'isfinite' not supported for the
                # input types, and the inputs could not be safely
                # coerced to any supported types according to the
                # casting rule ''safe''
                except TypeError:
                    pass

                last_best_energies = best_energies
                last_best_rmsds = best_rmsds
            mev += 1

    # run output
    update_all_fitness(ring)
    run_output(ring, f"last {mev}", save)
