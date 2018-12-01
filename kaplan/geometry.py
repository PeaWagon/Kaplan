
import vetee

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

def zmatrix_to_xyz(zmatrix):
    """Make xyz coordinates from a zmatrix com file."""
    # make a trivial com file with the zmatrix included
    # read in the com file with vetee.com(trivial_file)
    # make new vetee.xyz object
    # copy coords from zmatrix
    # return the object xyz coordinates:
    # list[[a1,x1,y1,z1], [a2,x2,y2,z2], ..., [an,xn,yn,zn]]
    # delete files
    xyz = []
    return xyz

def generate_zmatrix(parser, dihedrals):
    """Make a zmatrix comfile for a specific conformer.

    Parameters
    ----------
    parser : parser object
        Contains the original geometry for the molecule
        of interest.
    dihedrals : list(int)
        The dihedrals to be combined with
        the original geometry (in parser). 

    Returns
    -------
    zmatrix : str
        The full geometry specification of the
        molecule in a file of name zmatrix.

    """
    # generate a full geometry specification
    # in the form of a zmatrix using the dihedrals
    # found in the slot at pmem_index
    zmatrix = ""

    # generate zmatrix based on initial geom (parser)
    # replace the dihedrals with the dihedrals from the
    # pmem object
    # return a string
    return zmatrix
