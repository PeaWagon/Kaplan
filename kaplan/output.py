"""The purpose of this module is to handle writing output
after the conformer search has concluded.

This module decides where to put the output files, what
the output files will be called, and what the output
files contain.

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

Some new features that I'd like to add are
written as comments below."""

import os
import pickle

import numpy as np

from vetee.xyz import Xyz
from vetee.gaussian_options import periodic_table

from kaplan.geometry import get_geom_from_dihedrals
from kaplan.inputs import Inputs

# OUTPUT_FORMAT = 'xyz'

# FEATURES TODO:
# add option to change output format
# make some images representing the population using matplotlib
# stats file should include energies for each conformer and rmsd for each pair


def get_output_dir(loc="pwd"):
    """Determine the name of the output directory.

    Parameters
    ----------
    loc : str
        The parent directory to use as output.
        Defaults to "pwd", which means use the
        present working directory (current working
        directory). Another possible option is
        "home", which puts the output in the
        home directory (i.e. /user/home), but this
        option is only available for Linux users.
        If loc is not home or pwd, then the
        output will be generated in the given
        directory.

    Raises
    ------
    FileNotFoundError
        The user gave a location that does not exist.
    
    Notes
    -----
    The output is placed in kaplan_output under a job
    number, formatted as follows:
    loc/kaplan_output/job_0 # for the first job
    loc/kaplan_output/job_1 # for the second job
    etc.

    Returns
    -------
    output_dir : str
        The directory where the job output will
        be written.

    """
    if loc == "pwd":
        output_dir = os.path.join(os.getcwd(), "kaplan_output")
    elif loc == "home":
        output_dir = os.path.join(os.path.expanduser("~"), "kaplan_output")
    else:
        output_dir = os.path.join(os.path.abspath(loc), "kaplan_output")

    # first check that there is a place to put the
    # output files
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
        output_dir = os.path.join(output_dir, "job_0")
        os.mkdir(output_dir)
        return output_dir
    # iterate over existing jobs to determine
    # dir_num for newest job
    dir_contents = os.scandir(output_dir)
    # keep track of how many jobs have been run
    dir_nums = []
    for val in dir_contents:
        val = val.name.split("_")
        if val[0] == "job" and len(val) == 2:
            try:
                dir_nums.append(int(val[1]))
            except ValueError:
                pass
    new_dir = "job_" + str(max(dir_nums) + 1)
    output_dir = os.path.join(output_dir, new_dir)
    os.mkdir(output_dir)
    return output_dir


def run_output(ring, save, output_dir="pwd"):
    """Run the output module.

    Parameters
    ----------
    ring : object
       The final ring data structure after evolution.
    save : bool
        If True, writes a pickle file of the ring and
        input objects.
        If False, does not write any pickle files.
    output_dir : str
        The directory where kaplan_output is generated.
        Defaults to the current/present working directory.

    """
    inputs = Inputs()

    # find average fitness
    # find best pmem and its ring index
    total_fit = 0
    best_pmem = 0
    best_fit = 0
    for pmem in ring:
        if pmem is not None:
            fitg = pmem.fitness
            total_fit += fitg
            if fitg > best_fit:
                best_pmem = pmem.ring_loc
                best_fit = fitg
    average_fit = total_fit / ring.num_filled

    # generate and get output directory
    output_dir = get_output_dir()

    # write a stats file
    with open(os.path.join(output_dir, "stats-file.txt"), "w") as fout:
        fout.write(f"num conformers: {inputs.num_geoms}\n")
        fout.write(f"average fitness: {average_fit}\n")
        fout.write(f"best fitness: {best_fit}\n")
        fout.write(f"final percent filled: {100*ring.num_filled/ring.num_slots}%\n")
        fout.write("--------------------------------------------------------\n")
        for conf in range(inputs.num_geoms):
            fout.write(f"Conformer {conf}:\n")
            fout.write(f"Dihedrals: {ring[best_pmem].dihedrals[conf]}\n")
            fout.write(f"Energy: {ring[best_pmem].energies[conf]}\n")
            fout.write("--------------------------------------------------------\n")
        fout.write(f"RMSDs:\n")
        for rmsd in ring[best_pmem].rmsds:
            fout.write(f"Conf{rmsd[0]} & Conf{rmsd[1]} => {rmsd[2]}\n")

    # generate the output file for the best pmem
    for i, geom in enumerate(ring[best_pmem].all_coords):
        xyz = Xyz()
        # convert coords to angstroms
        for j, atom in enumerate(geom):
            # need to convert from numpy's int64 to regular python int
            label = periodic_table(int(inputs.atomic_nums[j]))
            xyz.coords.append([label, atom[0], atom[1], atom[2]])
        xyz.num_atoms = inputs.num_atoms
        xyz.comments = f"conformer {i}; energy {ring[best_pmem].energies[i]}"
        xyz.write_xyz(os.path.join(output_dir, f"conf{i}.xyz"))
    
    # pickel if requested
    if save:
        with open(os.path.join(output_dir,"ring.pickle"), "wb") as f:
            pickle.dump(ring, f)
        with open(os.path.join(output_dir, "inputs.pickle"), "wb") as f:
            pickle.dump(inputs, f)
