"""This module is responsible for dealing with converting
geometries from z-matrix format to xyz format. It relies
on the external library, Vetee, to do most of the conversions,
along with the python wrapper for openbabel, pybel."""

import vetee
import openbabel
import pybel

from kaplan.inputs import Inputs


class GeometryError(Exception):
    """Error raises when one of the geometry functions fails."""


def get_zmatrix_template(parser):
    """Make a zmatrix from the original geometry.

    Parameters
    ----------
    parser : object
        Parser object from vetee.

    Returns
    -------
    zmatrix : str
        The zmatrix (gzmat self.parserformat) for the parser
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
    if parser.comments is not None:
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
    mol = pybel.readstring("gzmat", zmatrix)
    for atom in mol.atoms:
        xyz.append([vetee.gaussian_options.periodic_table(atom.atomicnum),
                    atom.coords[0], atom.coords[1], atom.coords[2]])
    return xyz
