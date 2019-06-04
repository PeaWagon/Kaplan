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

from vetee.tools import periodic_table
from vetee.coordinates import write_xyz

from kaplan.inputs import Inputs

# OUTPUT_FORMAT = 'xyz'

# FEATURES TODO:
# add option to change output format
# make some images representing the population using matplotlib
# stats file should include energies for each conformer and rmsd for each pair





def run_output(ring, save):
    """Run the output module.

    Parameters
    ----------
    ring : object
       The final ring data structure after evolution.
    save : bool
        If True, writes a pickle file of the ring and
        input objects.
        If False, does not write any pickle files.

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
    output_dir = inputs.output_dir

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
        coords = []
        # convert coords to angstroms
        for j, atom in enumerate(geom):
            # need to convert from numpy's int64 to regular python int
            label = periodic_table(int(inputs.atomic_nums[j]))
            coords.append([label, atom[0], atom[1], atom[2]])
        comments = f"conformer {i}; energy {ring[best_pmem].energies[i]}"
        write_xyz({"_coords": coords, "_comments": comments},
                  os.path.join(output_dir, f"conf{i}.xyz"))
    
    # pickel if requested
    # note: cannot pickle obmol => TypeError: can't pickle SwigPyObject objects
    if save:
        with open(os.path.join(output_dir,"ring.pickle"), "wb") as f:
            pickle.dump(ring, f)
        with open(os.path.join(output_dir, "inputs.pickle"), "wb") as f:
            inputs.obmol = None
            pickle.dump(inputs, f)
