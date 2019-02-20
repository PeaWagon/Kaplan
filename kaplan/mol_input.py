"""
This module reads in the mol_input_file and
verifies the types of the input, the values
of the input, and the completeness of the
input.
"""

from kaplan.geometry import generate_parser
from kaplan.energy import BasisSetError, MethodError, run_energy_calc, DEFAULT_INPUTS

NUM_MOL_ARGS = 7
NUM_GA_ARGS = 12


def read_mol_input(mol_input_file):
    """Read in a mol input file.

    Parameters
    ----------
    mol_input_file : str
        The full path and name of the file to read.

    Raises
    ------
    FileNotFoundError
        Invalid path and file name.
    ValueError
        Incorrect number of arguments (or
        number of arguments is the same
        as the number of ga arguments -
        the order for the input files is
        switched).

    Returns
    -------
    mol_input_dict : dict
        The molecular input parameters organised
        as a dictionary.

    """
    mol_input_dict = dict()
    num_args = 0
    try:
        with open(mol_input_file, 'r') as fout:
            for line in fout:
                # remove \n char and separate keys and values
                # make everything lowercase to avoid
                # non-matching strings
                line = line[:-1].split(' = ')
                # ignore blank lines/long input
                if len(line) != 2 and line != "":
                    print(f"Warning: line - {line} - was ignored from the mol_input_file.")
                    continue
                mol_input_dict[line[0].lower()] = line[1]
                num_args += 1
                # go through each line and pull data and key
    except FileNotFoundError:
        raise FileNotFoundError("No such mol_input_file.")
    if num_args == NUM_GA_ARGS:
        raise ValueError("Wrong order of input files.")
    if num_args != NUM_MOL_ARGS:
        raise ValueError("Incorrect number of molecule arguments.")
    return mol_input_dict


def verify_mol_input(mol_input_dict):
    """Check molecular inputs are within bounds.

    Parameters
    ----------
    mol_input_dict : dict
        The information gathered from the input file.

    Returns
    -------
    parser : obj
        Vetee parser object.

    """
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

    # convert all but smiles string to lowercase
    for key in mol_input_dict:
        if key != "struct_input":
            # string is in case dictionary has
            # already been converted to integers/floats
            # pre-emptively
            mol_input_dict[key] = str(mol_input_dict[key]).lower()
    if mol_input_dict["struct_type"] != "smiles":
        mol_input_dict["struct_type"] = mol_input_dict["struct_type"].lower()

    # only program currently available is psi4
    # if a program is added, add default options in the energy module
    assert mol_input_dict["prog"] in DEFAULT_INPUTS

    # make sure the inputs are of the correct format
    assert mol_input_dict['struct_type'] in ('smiles', 'com', 'xyz', 'glog', "name", "cid")
    # check the structure file exists (if applicable)
    if mol_input_dict["struct_type"] in ("xyz", "glog", "com"):
        try:
            with open(mol_input_dict["struct_input"], "r"):
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
    except Exception:
        raise ValueError("Error when generating Parser object. Check the struct_input value.")

    # check initial geometry for convergence
    # check method and basis are in chosen program
    # if you add extra options, include them in the test here
    opt = {"prog": mol_input_dict["prog"],
           "qcm": mol_input_dict["qcm"],
           "basis": mol_input_dict["basis"],
           "charge": parser.charge,
           "multip": parser.multip
    }
    try:
        run_energy_calc(parser.coords, opt)
    except BasisSetError:
        raise BasisSetError(f"Invalid basis set: {mol_input_dict['basis']}")
    except QCMError:
        raise QCMError(f"Invalid quantum chemical method: {mol_input_dict['qcm']}")
    
    # if no error message, initial geometry converges, we are good
    return parser
