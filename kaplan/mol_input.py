
from kaplan.geometry import generate_parser
from kaplan.energy import prep_psi4_geom, run_energy_calc, check_psi4_inputs

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
    # make sure dict is non-empty
    assert len(mol_input_dict)
    # make sure it is actually a dictionary
    assert isinstance(mol_input_dict, dict)
    # key names
    key_names = {'qcm', 'basis', 'struct_input', 'struct_type', 'prog',
                 'charge', 'multip'}
    for key in key_names:
        try:
            assert key in mol_input_dict
        except AssertionError:
            raise ValueError(f"Misspelled/incorrectly formatted mol input parameter: {key}.")
    # ensure program used is psi4
    assert mol_input_dict["prog"] == "psi4"
    # check method and basis are in psi4
    try:
        check_psi4_inputs(mol_input_dict["qcm"], mol_input_dict["basis"])
    except ValueError:
        raise ValueError(f"Invalid basis set and/or method for psi4: {mol_input_dict['basis'], mol_input_dict['qcm']}")
    # make sure the inputs are of the correct format
    assert mol_input_dict['struct_type'] in ('smiles', 'com', 'xyz', 'glog', "name", "cid")
    # check the structure file exists (if applicable)
    if mol_input_dict["struct_type"] in ("xyz", "glog", "com"):
        try:
            with open(mol_input_dict["struct_input"], 'r') as f:
                pass
        except FileNotFoundError:
            raise FileNotFoundError("No such struct_input file.")
    # check the charge and multiplicity are integers
    try:
        mol_input_dict["multip"] = int(mol_input_dict["multip"])
        mol_input_dict["charge"] = int(mol_input_dict["charge"])
    except ValueError:
        raise ValueError(f"Charge and multiplicity should be integer values.")
    assert mol_input_dict["multip"] > 0
    # try to make the struct object using vetee
    try:
        parser = generate_parser(mol_input_dict)
    except Exception as e:
        raise ValueError("Error when generating Parser object. Check the struct_input value.")
    # check here if error message is raised
    geom = prep_psi4_geom(parser.coords, parser.charge, parser.multip)
    run_energy_calc(geom, mol_input_dict["qcm"], mol_input_dict["basis"])
    # if no error message, initial geometry converges, we are good
    return parser
    
