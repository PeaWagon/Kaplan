
from kaplan.geometry import generate_parser
from kaplan.energy import prep_psi4_geom, run_energy_calc

"""
This module reads in the mol_input_file and
verifies the types of the input, the values
of the input, and the completeness of the
input.
"""

NUM_MOL_ARGS = 7

def read_mol_input(mol_input_file):
    mol_input_dict = dict()
    num_args = 0
    try:
        with open(mol_input_file, 'r') as f:
            for line in f:
                # remove \n char and separate keys and values
                # make everything lowercase to avoid
                # non-matching strings
                line = line[:-1].lower().split(' = ')
                # ignore blank lines/long input
                if len(line) != 2:
                    print(f"Warning: line - {line} - was ignored from the mol_input_file.")
                    continue
                mol_input_dict[line[0]] = line[1]
                num_args += 1
                # go through each line and pull data and key
    except FileNotFoundError:
        raise FileNotFoundError("No such mol_input_file.")
    if num_args != NUM_MOL_ARGS:
        raise ValueError("Incorrect number of molecule arguments.")
    return mol_input_dict

def verify_mol_input(mol_input_dict):
    assert len(mol_input_dict)
    assert isinstance(mol_input_dict, dict)
    # make sure the inputs are of the correct format
    assert mol_input_dict['struct_type'] in ('smiles', 'com', 'xyz', 'glog')
    # try to make the struct object using vetee
    # somehow run an energy calculation to test its convergence
    # and to ensure that the method and basis set are available
    # in the selected program
    parser = generate_parser(mol_input_dict)
    # check here if error message is raised
    run_energy_calc(prep_psi4_geom(parser))
    # if no error message, initial geometry converges, we are good
    return parser
    
