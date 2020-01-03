"""The purpose of this module is to handle writing output
after the conformer search has concluded, and to allow
the user to keep track of the state of the population
during evolution.

The inputs module decides where to put the output files.
This module decides what the output files will be called,
and what the output files contain.

The default output directory is:
os.cwd/kaplan_output/

Each time Kaplan is run, a new job number is created:
os.cwd/kaplan_output/job_0
os.cwd/kaplan_output/job_1
        ....
os.cwd/kaplan_output/job_n-1
where n is the number of times Kaplan is run.

Note: if run directories are deleted, Kaplan
will generate new job numbers depending on the
highest value it finds (so if job 0, 1, 5 and 10
are present, then the next job will be job_11).

"""

import os
import pickle
from time import asctime

from kaplan.inputs import Inputs, InputError
from kaplan.geometry import update_obmol, write_coords, periodic_table

import kaplan.tools as kt

import pkg_resources


# enables use of kaplan.__version__
__version__ = pkg_resources.get_distribution("kaplan").version


def run_output(ring, current_mev, save):
    """Run the output module.

    Parameters
    ----------
    ring : object
       The current state of the ring data structure.
       Assumes that the fitness values are up to date.
    current_mev : int, str
        Current mating event. If the current mev is
        smaller than the number of mating events total,
        then only the stats file and the ring object
        should be updated. If set to "last", then
        the inputs object is written (minus obmol
        object, as swig objects cannot be pickled).
    save : bool
        If True, writes a pickle file of the ring and
        input objects.
        If False, does not write any pickle files.

    Notes
    -----
    If the final output is run, this will delete the
    obmol object from the inputs. Therefore, do not
    call the final output unless there is nothing
    more to be done with obmol.

    """
    # make sure no other keywords are used for mev
    assert isinstance(current_mev, int) or "last" in current_mev

    inputs = Inputs()

    # find best pmem and its ring index
    best_pmem, best_fit = ring.best_pmem
    if best_pmem is None:
        raise InputError("Ring does not have a best pmem. Try updating fitness values.")

    units = kt.units_by_prog[inputs.prog]

    # do the best pmem files before the stats file update
    # so that the final energies and RMSD values are reported
    # based on a major optimisation
    # generate the output file for the best pmem
    if not isinstance(current_mev, int):
        major_coords = ring[best_pmem].setup(major=True)
        for i, coords in enumerate(major_coords):
            if coords is None:
                print("Warning: could not generate an output file.")
                print(f"No valid coordinates for conformer {i}.")
            else:
                outfile = os.path.join(inputs.output_dir, f"conf{i}.xyz")
                comments = f"{inputs.name} | "
                comments += f"conformer {i} | energy {ring[best_pmem].energies[i]}"
                comments += f" ({units})"
                write_coords(coords, inputs.atomic_nums, outfile, comments)
                # make a 3d image of each conformer
                kt.plot_3d(outfile)
        # make energy plot for best pmem
        try:
            kt.energy_barplot(ring[best_pmem], inunits=units)
        # ValueError: Image size of 73500x1800 pixels is too large.
        # It must be less than 2^16 in each direction.
        except ValueError:
            print("Warning: too many conformers in pmem.")
            print("Figure exceeded max allowed size.")
            print("Cannot make output plot.")
        # if more than one geom per pmem, make energy delta/rmsd plot
        if inputs.num_geoms > 1:
            try:
                kt.energy_rmsd_scatter(ring[best_pmem], inunits=units)
            # pmem has no energies
            except InputError:
                print("Warning. Best pmem has no valid energy values.")
                print("Cannot make output plot.")
            # ValueError: Image size of 73500x1800 pixels is too large.
            # It must be less than 2^16 in each direction.
            except ValueError:
                print("Warning: too many conformers in pmem.")
                print("Figure exceeded max allowed size.")
                print("Cannot make output plot.")
        # if at least 2 dihedrals, make a dihedral heatmap
        if inputs.num_diheds > 1:
            kt.dihedrals_heatmap(ring)

    # write a stats file
    with open(os.path.join(inputs.output_dir, "stats-file.txt"), "a") as fout:
        fout.write("--------------------------------------------------------\n")
        mev = f"MEV {current_mev}"
        fout.write(f"{mev:^55}\n")
        fout.write(f"{asctime():^55}\n")
        kv = f"Kaplan Version: {__version__}"
        fout.write(f"{kv:^55}\n")
        fout.write("--------------------------------------------------------\n")
        fout.write(f"average fitness:         {ring.mean_fitness} +/- {ring.stdev_fitness}\n")
        fout.write(f"median fitness:          {ring.median_fitness}\n")
        fout.write(f"best fitness:            {best_fit} (pmem{best_pmem})\n")
        fout.write(f"num conformers per pmem: {inputs.num_geoms}\n")
        fout.write(f"slots filled:            {ring.num_filled}/{ring.num_slots}\n")
        fout.write(f"average energy ({units}): {ring.mean_energy} +/- {ring.stdev_energy}\n")
        fout.write(f"median energy ({units}):  {ring.median_energy}\n")
        if inputs.num_geoms > 1:
            fout.write(f"average rmsd:            {ring.mean_rmsd} +/- {ring.stdev_rmsd}\n")
            fout.write(f"median rmsd:             {ring.median_rmsd}\n")
        fout.write("--------------------------------------------------------\n")
        for conf in range(inputs.num_geoms):
            fout.write(f"Conformer: {conf}\n")
            fout.write(f"Dihedrals: {ring[best_pmem].dihedrals[conf]}\n")
            fout.write(f"Energy ({units}): {ring[best_pmem].energies[conf]}\n")
            fout.write("--------------------------------------------------------\n")
        fout.write(f"RMSDs (excluding {inputs.exclude_from_rmsd}):\n")
        for rmsd in ring[best_pmem].rmsds:
            fout.write(f"Conf{rmsd[0]} & Conf{rmsd[1]} => {rmsd[2]}\n")
        fout.write("\n")

    # pickel if requested
    # note: cannot pickle obmol => TypeError: can't pickle SwigPyObject objects
    if save:
        with open(os.path.join(inputs.output_dir, "ring.pickle"), "wb") as f:
            pickle.dump(ring, f)
        if not isinstance(current_mev, int):
            with open(os.path.join(inputs.output_dir, "inputs.pickle"), "wb") as f:
                inputs.obmol = None
                pickle.dump(inputs, f)
