
import vetee
import openbabel
import pybel

class GeometryError(Exception):
    """Error raises when one of the geometry functions fails."""
    pass

def generate_parser(mol_input_dict):
    """Returns parser object (from vetee).

    Parameters
    ----------
    mol_input_dict : dict
        Need following keys:
        1. struct_type
        2. struct_input
        3. qcm
        4. basis
        5. charge
        6. multip

    """
    try:
        if mol_input_dict['struct_type'] == 'xyz':
            parser = vetee.xyz.Xyz(mol_input_dict['struct_input'])
        elif mol_input_dict['struct_type'] == 'com':
            parser = vetee.com.Com(mol_input_dict['struct_input'])
        elif mol_input_dict['struct_type'] == 'glog':
            parser = vetee.glog.Glog(mol_input_dict['struct_input'])
        elif mol_input_dict['struct_type'] in ('smiles', 'cid', 'name'):
            parser = vetee.structure.Structure(mol_input_dict['struct_type'], mol_input_dict['struct_input'])
        else:
            raise GeometryError(f"The struct_type input {mol_input_dict['struct_type']} is not available.")
    except Exception as e:
        print(e)
        raise GeometryError(f"Unable to make a parser object. {mol_input_dict['struct_type']}: {mol_input_dict['struct_input']}")
    # update basis set and method
    # note that the mol_input_file has precedence over
    # the contents of struct_input
    parser.parse_gkeywords(f"#{mol_input_dict['qcm']} {mol_input_dict['basis']}")
    # update charge and multiplicity, again mol_input_file
    # takes precedence over whatever was decided when the
    # file was read by vetee
    if not isinstance(mol_input_dict['charge'], int) or not isinstance(mol_input_dict['multip'], int):
        raise GeometryError("Charge and multiplicity should be integers. Unable to make a parser object.")
    parser.charge = mol_input_dict['charge']
    parser.multip = mol_input_dict['multip']
    return parser

def get_zmatrix_template(parser):
    """Make a zmatrix from the original geometry.

    Parameters
    ----------
    parser : object
        Parser object from vetee.

    Returns
    -------
    zmatrix : str
        The zmatrix (gzmat format) for the parser
        molecule. Will be used to combine with
        dihedral angles.

    """
    # from vetee
    obmol = openbabel.OBMol()
    # add coordinates for each atom
    for atom in parser.coords:
        obatom = openbabel.OBAtom()
        atomicnum = vetee.gaussian_options.periodic_table(atom[0])
        obatom.SetAtomicNum(atomicnum)
        obatom.SetVector(atom[1], atom[2], atom[3])
        obmol.AddAtom(obatom)
    # set charge, multiplicity, and comments (title)
    obmol.SetTotalCharge(parser.charge)
    obmol.SetTotalSpinMultiplicity(parser.multip)
    obmol.SetTitle(parser.comments)
    # convert the obmol to a pybel Molecule
    pybelmol = pybel.Molecule(obmol)
    # generate a zmatrix string using the obmol input
    zmatrix = pybelmol.write("gzmat")
    return zmatrix

def update_zmatrix(zmatrix, dihedrals):
    """Make a new zmatrix with given dihedral angles.

    Parameters
    ----------
    zmatrix : str
        Contains the original geometry for the molecule
        of interest.
    dihedrals : list(int)
        The dihedrals to be combined with
        the original geometry. 

    Returns
    -------
    zmatrix : str
        The full geometry specification of the
        molecule in a zmatrix.

    """
    zmatrix_list = zmatrix.split('\n')
    dihedral_num = 0
    for i, line in enumerate(zmatrix_list):
        if line.startswith('d') and '=' in line:
            line = line[:line.index('=')+1] + str(dihedrals[dihedral_num])
            zmatrix_list[i] = line
            dihedral_num += 1
    new_zmatrix = '\n'.join(zmatrix_list)
    return new_zmatrix

def zmatrix_to_xyz(zmatrix):
    """Make xyz coordinates from a zmatrix string.

    Parameters
    ----------
    zmatrix : str
        A zmatrix for the input molecule.

    Returns
    -------
    xyz : list(list(str, float, float, float))
        The xyz coordinates of the molecule
        as converted from the zmatrix, including
        atom types.
    [[a1,x1,y1,z1], [a2,x2,y2,z2], ..., [an,xn,yn,zn]]

    """
    xyz = []
    mol = pybel.readstring("gzmat", new_zmatrix)
    for atom in mol.atoms:
        xyz.append([vetee.gaussian_options.periodic_table(atom.atomicnum),
                    atom.coords[0], atom.coords[1], atom.coords[2]])
    return xyz

