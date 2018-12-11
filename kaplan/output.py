
import os

from vetee.xyz import Xyz

from kaplan.geometry import update_zmatrix, zmatrix_to_xyz

# OUTPUT_FORMAT = 'xyz'

# directory for this test file
file_dir = os.path.dirname(os.path.realpath(__file__))
OUTPUT_DIR = os.path.join(file_dir, "kaplan_output")

# to change the output location, change this variable to
# the desired directory and uncomment the line
# note: the output may not work if multiple subdirectories
# need to be generated
#OUTPUT_DIR = "~/kaplan_output"

def run_output(ring):
    """Run the output module.

    Parameters
    ----------
    ring : object
       the final ring data structure after evolution.

    """
    # first check that there is a place to put the output files
    if not os.path.isdir(OUTPUT_DIR):
        os.mkdir(OUTPUT_DIR)

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

    # write a stats file
    with open(os.path.join(OUTPUT_DIR, "stats-file.txt"), "a") as f:
        f.write(f"+++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
        f.write(f"average fitness: {average_fit}\n")
        f.write(f"best fitness: {best_fit}\n")
        f.write(f"final percent filled: {100*ring.num_filled/ring.num_slots}%\n")

    # generate the output file for the best pmem
    # TODO: figure out a way to make the files specific to the input
    # molecule (right now they overwrite each other)
    # write the xyz files
    for geom in range(ring.num_geoms):
        xyz_coords = zmatrix_to_xyz(update_zmatrix(ring.zmatrix, ring[best_pmem].dihedrals[geom]))
        xyz = Xyz()
        xyz.coords = xyz_coords
        xyz.num_atoms = ring.num_atoms
        xyz.comments = f"conformer {geom}"
        xyz.write_xyz(os.path.join(OUTPUT_DIR, f"conf{geom}.xyz"))
        


