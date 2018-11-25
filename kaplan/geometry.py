

from kaplan.energy import prep_psi4_geom, run_energy_calc


def generate_parser(mol_input_dict):
    """Returns parser object (from vetee)."""
    if mol_input_dict['struct_type'] == 'xyz':
        parser = vetee.xyz.Xyz(mol_input_dict['struct_input'])
    elif mol_input_dict['struct_type'] == 'com':
        parser = vetee.com.Com(mol_input_dict['struct_input'])
    elif mol_input_dict['struct_type'] == 'glog':
        parser = vetee.glog.Glog(mol_input_dict['struct_input'])
    elif mol_input_dict['struct_type'] in ('smiles', 'cid', 'name'):
        parser = vetee.structure.Structure(mol_input_dict['struct_type'], mol_input_dict['struct_input'])
    else:
        raise NotImplementedError(f"The struct_type input {mol_input_dict['struct_type']} is not yet supported.")
    # update basis set and method
    # note that the mol_input_file has precedence over
    # the contents of struct_input
    parser.parse_gkeywords(f"#{mol_input_dict['qcm']} {mol_input_dict['basis']}")
    # update charge and multiplicity, again mol_input_file
    # takes precedence over whatever was decided when the
    # file was read by vetee
    parser.charge = mol_input_dict['charge']
    parser.multip = mol_input_dict['multip']
    return parser

def write_output(dihedrals):
    pass

def generate_zmatrix(parser, new_dihedrals):
    """Make a zmatrix string.

    Parameters
    ----------
    parser : parser object
        Contains the original geometry for the molecule
        of interest.
    new_dihedrals : list of floats
        The new dihedral angles for the molecule.

    """
    pass
