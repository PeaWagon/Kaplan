"""The purpose of this module is to handle writing output
after the conformer search has concluded.

This module decides where to put the output files, what
the output files will be called, and what the output
files contain.

The output directory parent is:
os.cwd/kaplan_output/

Each time Kaplan is run, a new job number is created:
os.cwd/kaplan_output/job_0
os.cwd/kaplan_output/job_1
        ....
os.cwd/kaplan_output/job_n
where n is the number of times Kaplan is run.

Note: if run directories are deleted, Kaplan
will generate new job numbers depending on the
highest value it finds (so if job 0, 1, 5 and 10
are present, then the next job will be job_11).

Some new features that I'd like to add are
written as comments below."""

import os

from vetee.xyz import Xyz

from kaplan.geometry import update_zmatrix, zmatrix_to_xyz

# OUTPUT_FORMAT = 'xyz'

# FEATURES TODO:
# add option to change output format
# add option to write to a user-specified output directory
# make some images representing the population using matplotlib


def get_output_dir(loc="pwd"):
    """Determine the name of the output directory.

    Parameters
    ----------
    loc : str
        The parent directory to use as output.
        Defaults to "pwd", which means use the
        present working directory (current working
        directory). The other possible option is
        "home", which puts the output in the
        home directory (i.e. /user/home), but this
        option is only available for Linux users.

    Raises
    ------
    AssertionError
        The user gave an option that was not "home"
        or "pwd".

    Returns
    -------
    output_dir : str
        The directory where the job output will
        be written.

    """
    assert loc in ["pwd", "home"]
    output_dir = os.path.join(os.getcwd(),
                              "kaplan_output")
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


def run_output(ring):
    """Run the output module.

    Parameters
    ----------
    ring : object
       The final ring data structure after evolution.

    """
    # find average fitness
    # find best pmem and its ring index
    total_fit = 0
    best_pmem = 0
    best_fit = 0
    for pmem in ring.pmems:
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
        fout.write(f"num conformers: {ring.num_geoms}\n")
        fout.write(f"average fitness: {average_fit}\n")
        fout.write(f"best fitness: {best_fit}\n")
        fout.write(f"final percent filled: {100*ring.num_filled/ring.num_slots}%\n")

    # generate the output file for the best pmem
    for geom in range(ring.num_geoms):
        xyz_coords = zmatrix_to_xyz(update_zmatrix(ring.zmatrix, ring[best_pmem].dihedrals[geom]))
        xyz = Xyz()
        xyz.coords = xyz_coords
        xyz.num_atoms = ring.num_atoms
        xyz.comments = f"conformer {geom}"
        xyz.write_xyz(os.path.join(output_dir, f"conf{geom}.xyz"))
