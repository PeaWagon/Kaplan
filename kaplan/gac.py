
from kaplan.ga_input import read_ga_input, verify_ga_input

"""

In charge of evolving a population for the set
number of mating events.

Can also decide upon another stopping condition
whereby there is a set of conformers with a suitably
low energy (i.e. good fitness).

Haven't decided on whether to make this an object or
not yet.

Needs to contain:
* stopping condition
* what kinds of mutations to perform
* selection rules (i.e. tournament)

"""

def run_kaplan(ga_input_file, mol_input):
    """Run the Kaplan programme.

    Parameters
    ----------
    ga_input_file : str
        The input file containing genetic algorithm constants.
    mol_input : object (vetee)
        A vetee Parser object (or child class).

    """
    ga_input_dict = read_ga_input(ga_input_file)
    verify_ga_input(ga_input_dict)
