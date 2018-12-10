
import sys

from kaplan.ga_input import read_ga_input, verify_ga_input
from kaplan.mol_input import read_mol_input, verify_mol_input
from kaplan.ring import Ring, RingEmptyError
from kaplan.tournament import run_tournament
from kaplan.output import run_output

"""

In charge of evolving a population for the set
number of mating events.

"""

def run_kaplan(ga_input_file, mol_input_file):
    """Run the Kaplan programme.

    Parameters
    ----------
    ga_input_file : str
        The input file containing genetic algorithm constants.
    mol_input_file : str
        The input file containing the molecular information.

    """
    # read in and verify ga_input_file
    ga_input_dict = read_ga_input(ga_input_file)
    verify_ga_input(ga_input_dict)

    # read in and verify mol_input_file
    # check that initial geometry converges
    # and construct a parser object
    mol_input_dict = read_mol_input(mol_input_file)
    parser = verify_mol_input(mol_input_dict)

    # check that inputs agree on a very trivial level
    assert ga_input_dict['num_atoms'] == len(parser.coords)

    # make a ring
    ring = Ring(ga_input_dict['num_geoms'], 
                ga_input_dict['num_atoms'],
                ga_input_dict['num_slots'],
                ga_input_dict['pmem_dist'],
                ga_input_dict['fit_form'],
                ga_input_dict['coef_energy'],
                ga_input_dict['coef_rmsd'],
                parser)

    # fill ring with an initial population
    ring.fill(ga_input_dict['num_filled'], 0)

    # run the mevs
    for mev in range(ga_input_dict['num_mevs']):
        try:
            print(mev)
            run_tournament(ga_input_dict['t_size'],
                           ga_input_dict['num_muts'],
                           ga_input_dict['num_swaps'],
                           ring, mev)
        except RingEmptyError:
            ring.fill(ga_input_dict['num_filled'], mev)

    # run output
    run_output(ring)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise FileNotFoundError("Please include the ga_input_file and the mol_input_file as arguments to the program.")
    run_kaplan(sys.argv[1], sys.argv[2])
