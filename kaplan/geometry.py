"""This module is responsible for dealing with converting
geometries from z-matrix format to xyz format. It relies
on the external library, Vetee, to do most of the conversions,
along with the python wrapper for openbabel, pybel.

New changes to how dihedrals are used are being made.

Now this program will use saddle (GOpt) to generate
dihedral angles.

"""

import vetee
import openbabel
import pybel
import numpy as np
import os

from saddle.internal import Internal
from saddle.coordinate_types import DihedralAngle

#from kaplan.inputs import Inputs

# 1 "angstrom" = 1.8897261339213 "atomic unit (length)"
#AU = 1.8897261339213
# https://github.com/psi4/psi4/blob/fbb2ff444490bf6b43cb6e027637d8fd857adcee/psi4/include/psi4/physconst.h
# convert angstroms to atomic units (or bohr from psi4)
ANG_TO_AU = lambda x : x/0.52917721067
# values for dihedral angles in radians
MIN_VALUE = -2*np.pi
MAX_VALUE = 2*np.pi
# convert radians to degrees for __str__ method
RAD_TO_DEGREES = lambda x : x*180/np.pi


def get_struct_info(struct_input, struct_type="name", prog="psi4"):
    """Get structure information from an input query using Vetee.

    struct_input : str
        The value of the input to be parsed. Examples
        include files (with full path), CID value, and
        molecule name.
    struct_type : str
        One of name, cid, smiles, inchi, inchikey, com,
        or xyz. Specifies type of input provided. Defaults
        to name.
    prog : str
        Setup structural information according to program.
    
    Returns
    -------
    struct_info : vetee.job.Job object instance
        Specifies all information needed to write a molecule.
    
    """
    AVAIL_STRUCTS = (
        "xyz", "com", "smiles", "inchi", "inchikey",
        "name", "cid",
    )
    struct_type = struct_type.lower()
    assert struct_type in AVAIL_STRUCTS
    if struct_type in ("com", "xyz"):
        assert os.path.isfile(struct_input)
    
    struct_info = vetee.job.Job(struct_type, struct_input, prog)
    try:
        # kaplan=True to avoid not-implemented error
        struct_info.setup(kaplan=True)
    # Note: JobError raised when unable to set charge and multiplicity
    except vetee.job.JobError:
        pass
    return struct_info


def create_obmol(xyz_file, charge, multip):
    """Use openbabel to create an OBMol object from xyz file.
    
    Parameters
    ----------
    xyz_file : str
        Full file and path with which to generate OBMol.
        Units for the xyz file are angstroms.
    charge : int
        The overall charge of the molecule.
    multip : int
        The overall multiplicity of the molecule.

    Raises
    ------
    AssertionError
        xyz_file does not exist or is not a file.
    
    Returns
    -------
    obmol : openbabel.OBMol object
        Containing information from xyz file.
    
    """
    assert os.path.isfile(xyz_file)
    mol = pybel.readfile("xyz", xyz_file).__next__()
    mol = mol.OBMol
    if mol.GetTotalCharge() != charge:
        print("Warning: openbabel molecule charge different from input charge.")
        mol.SetTotalCharge(charge)
    if mol.GetTotalSpinMultiplicity() != multip:
        print("Warning: openbabel molecule multiplicity different from input multiplicity.")
        mol.SetTotalSpinMultiplicity(multip)
    return mol


def update_obmol(obmol, dihedrals, new_dihed):
    """Update an OBMol's torison angles.
    
    Parameters
    ----------
    obmol : openbabel.OBMol object
        The molecule to apply the changes to.
    dihedrals : list(tuple(int,int,int,int))
        List of tuples representing 4 atom
        indices for each dihedral.
    new_dihed : list(float) or np.array(num_dihed)
        New dihedrals to apply in radians.
        Should match index in dihedrals.
    
    Returns
    -------
    coords : np.array((num_atoms,3), float)
        New coordinates in angstroms.

    """
    assert len(dihedrals) == len(new_dihed)
    # obmol.BeginModify()
    for i, newt in enumerate(new_dihed):
        # get atom uses Idx which is atom index + 1
        atom1 = obmol.GetAtom(dihedrals[i][0]+1)
        atom2 = obmol.GetAtom(dihedrals[i][1]+1)
        atom3 = obmol.GetAtom(dihedrals[i][2]+1)
        atom4 = obmol.GetAtom(dihedrals[i][3]+1)
        # set torsion requires atom object instances
        # and input is in radians not degrees
        obmol.SetTorsion(atom1, atom2, atom3, atom4, newt)
    # obmol.EndModify()
    new_coords = np.zeros((obmol.NumAtoms(), 3), float)
    for i, atom in enumerate(openbabel.OBMolAtomIter(obmol)):
        new_coords[i][0] = atom.x()
        new_coords[i][1] = atom.y()
        new_coords[i][2] = atom.z()
    return new_coords


def set_coords(obmol, coords):
    """Set the coordinates for an openbabel molecule.

    Notes
    -----
    This function was written such that, if changes
    to dihedral angles broke atomic bonds, the original
    geometry could be maintained.

    """
    obmol.BeginModify()
    for i, atom in enumerate(openbabel.OBMolAtomIter(obmol)):
        # needs to be explicit (doesn't accept tuple, list, np.array)
        atom.SetVector(coords[i][0], coords[i][1], coords[i][2])
    obmol.EndModify()
    return obmol


def get_rings(obmol):
    """Get ring information from an OBMol object.

    Parameters
    ----------
    obmol : openbabel.OBMol object
        The molecule to search for rings.
    
    Returns
    -------
    rings : list(tuple(int, bool, str, list(int)))
        Each element in the list is information about
        a ring in the input molecule, stored as a
        four element tuple. The tuple describes:
        [0] how many atoms are in the ring?
        [1] is the ring aromatic?
        [2] what is the ring name (ex: benzene)?
        [3] what atom indices are involved in the ring?
    Note: [3] is connected in order around the ring (from _path).

    """
    rings = []
    for i, r in enumerate(openbabel.OBMolRingIter(obmol)):
        atoms = [x-1 for x in r._path] # convert Idx to index
        rings.append((r.Size(), r.IsAromatic(), r.GetType(), atoms))
    return rings


def get_atomic_nums(obmol):
    """Get list of atomic numbers for molecule."""
    atomic_nums = []
    for atom in openbabel.OBMolAtomIter(obmol):
        atomic_nums.append(atom.GetAtomicNum())
    return atomic_nums


def construct_fours(ring):
    """Construct all connected pairs of four atoms around a ring.

    Parameters
    ----------
    ring : tuple(int, bool, str, list(int))
        A ring contains a:
        [0] size (num atoms in ring).
        [1] True = is aromatic, False = not aromatic
        [2] name of ring type (ex: benzene)
        [3] atom indices in connected order around ring
    
    Raises
    ------
    AssertionError
        The ring is smaller than 4 items.

    Returns
    -------
    fours : list(tuple(int,int,int,int))
        All sets of 4 connected atoms in the ring (incl.
        clockwise and counter-clockwise directions).
    
    """
    ring_size = ring[0]
    assert ring_size > 3
    loop = ring[3] + ring[3][:3]
    fours = []
    for i in range(ring_size):
        piece = loop[i:i+4]
        piece2 = piece[::-1]
        fours.append(tuple(piece))
        fours.append(tuple(piece2))
    return fours


def remove_ring_dihed(rings, diheds):
    """Remove dihedral angles that are only in a ring.
    
    Parameters
    ----------
    rings : list(int, bool, str, list(int))
        As generated from get_rings function in this module.
    diheds : list(tuple(int,int,int,int))
        Each element is a tuple representing the 4 atom
        indices of the dihedral angle.
    
    Returns
    -------
    diheds : list(tuple(int,int,int,int))
        Updated set of dihedral angles, except with
        the angles that encompass only atoms in a ring
        removed.
    
    """
    to_remove = []
    for ring in rings:
        if ring[0] < 4:
            continue
        fours = construct_fours(ring)
        for dihed in diheds:
            if dihed in fours:
                to_remove.append(dihed)
    for remove in to_remove:
        diheds.remove(remove)
    return diheds




def get_min_dihed(xyz_file, charge, multip):
    """Get minimal dihedrals from file using GOpt.
    
    Parameters
    ----------
    xyz_file : str
        Full file and path with which to generate GOpt mol.
        Units for the xyz file are angstroms.
    charge : int
        Charge of the molecule.
    multip : int
        Multiplicity of the molecule.

    Raises
    ------
    AssertionError
        xyz_file does not exist or is not a file.
    
    Returns
    -------
    min_dihed : list(tuple(int,int,int,int))
        Tuple containing a tuple of 4 atom indices. The
        dihedral a-b-c-d as a tuple is (a,b,c,d).
    
    """
    assert os.path.isfile(xyz_file)
    mol = Internal.from_file(xyz_file, charge, multip)
    mol.auto_select_ic(minimum=True)
    min_diheds = []
    for intern_coord in mol.ic:
        if isinstance(intern_coord, DihedralAngle):
            atoms = intern_coord.atoms
            # GOpt atoms are np.int64 (not sure openbabel will like that)
            min_diheds.append(
                 (int(atoms[0]), int(atoms[1]),
                  int(atoms[2]), int(atoms[3])))
            #     float(intern_coord.value)) # value for dihedral angle
            
    return min_diheds
