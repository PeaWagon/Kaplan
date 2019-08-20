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

from copy import deepcopy

from saddle.internal import Internal
from saddle.coordinate_types import DihedralAngle
from saddle.errors import NotConvergeError


# values for dihedral angles in radians
# dihedrals should be positive values only
# so that mean can be calculated without
# having to change the dihedrals
# i.e. 3pi/4 and -pi/4 average should be
# 10pi/4 not pi/4 since -pi/4 is really
# 7pi/4
MIN_VALUE = 0
MAX_VALUE = 2 * np.pi

# 1 "angstrom" = 1.8897261339213 "atomic unit (length)"
# AU = 1.8897261339213
# https://github.com/psi4/psi4/blob/fbb2ff444490bf6b43cb6e027637d8fd857adcee/psi4/include/psi4/physconst.h
# convert angstroms to atomic units (or bohr from psi4)
geometry_units = {
    "angstroms": {
        "atomic units": lambda x: x / 0.52917721067,
        "angstroms": lambda x: x,
    },
    "atomic units": {
        "angstroms": lambda x: x * 0.52917721067,
        "atomic_units": lambda x: x,
    },
    "radians": {
        "degrees": lambda x: x * 180 / np.pi,
        "radians": lambda x: x,
    },
    "degrees": {
        "radians": lambda x: x * np.pi / 180,
        "degrees": lambda x: x,
    },
}


class GeometryError(Exception):
    """Raised when geometry problem occurs."""


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
    if struct_type in ("com", "xyz") and not os.path.isfile(struct_input):
        raise FileNotFoundError(f"No such input file: {struct_input}")
    struct_info = vetee.job.Job(struct_type, struct_input, prog)
    try:
        # kaplan=True to avoid not-implemented error
        struct_info.setup(kaplan=True)
    # Note: JobError raised when unable to set charge and multiplicity
    except vetee.job.JobError:
        pass
    return struct_info


def create_obmol_from_string(str_type, value):
    """Use openbabel to create an OBMol object from a string.

    Parameters
    ----------
    str_type : str
        Type of string being passed to the function.
        Should be one of: "smiles", "inchi".
    value : str
        The string with which to generate
        a molecule.

    Examples
    --------
    This function was written so that the user
    could create structures using SMILES strings
    that are not in Pubchem. Cyclooctatetraene (COT)
    was given a 2- charge using this method via
    the following smiles string:
        C1=C[CH-]=CC=C[CH-]=C1

    Notes
    -----
    Openbabel cannot construct molecules from
    InchiKey (example: GSYKYIXVPVFJRF-UHFFFAOYSA-N).

    Returns
    -------
    openbabel.OBMol object

    """
    if str_type == "smiles":
        str_type = "smi"
    try:
        mol = pybel.readstring(str_type, value)
    except OSError:
        raise GeometryError(f"Could not read {str_type} string with Openbabel: {value}") from None
    except ValueError:
        raise GeometryError(f"Openbabel does not support the input format: {str_type}") from None
    mol.make3D()
    return mol.OBMol


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
        print(f"Warning: openbabel molecule charge ({mol.GetTotalCharge()})\
                \ndifferent from input charge ({charge}).")
        mol.SetTotalCharge(charge)
    if mol.GetTotalSpinMultiplicity() != multip:
        print(f"Warning: openbabel molecule multiplicity ({mol.GetTotalSpinMultiplicity()})\
                \ndifferent from input multiplicity ({multip}).")
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
    new_dihed : list(float) or np.array(num_diheds)
        New dihedrals to apply in radians.
        Should match index in dihedrals.

    Returns
    -------
    coords : np.array((num_atoms,3), float)
        New coordinates in angstroms.

    """
    if len(dihedrals) != len(new_dihed):
        raise GeometryError(
            f"\nCurrent number of dihedrals ({len(dihedrals)})\
              \ndoes not match given number ({len(new_dihed)})"
        )
    # obmol.BeginModify()
    for i, newt in enumerate(new_dihed):
        # get atom uses Idx which is atom index + 1
        atom1 = obmol.GetAtom(dihedrals[i][0] + 1)
        atom2 = obmol.GetAtom(dihedrals[i][1] + 1)
        atom3 = obmol.GetAtom(dihedrals[i][2] + 1)
        atom4 = obmol.GetAtom(dihedrals[i][3] + 1)
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


def get_torsion(obmol, dihedral):
    """Get the torsion angle of an openbabel molecule.

    Parameters
    ----------
    obmol : openbabel.OBMol object
        The molecule for which to calculate the torsion
        angle.
    dihedral : tuple(int, int, int, int)
        The dihedral angle for which to calculate
        the angle. Each integer is an atom index.

    Notes
    -----
    By default, Openbabel's GetTorsion accepts IdX
    values, which means that the atom indices are +1
    from the input dihedral values. Also by default,
    the value of the dihedral is returned in degrees.

    Returns
    -------
    torsion angle in radians

    """
    return geometry_units["degrees"]["radians"](
        obmol.GetTorsion(dihedral[0] + 1,
                         dihedral[1] + 1,
                         dihedral[2] + 1,
                         dihedral[3] + 1)
    )


def get_coords(obmol):
    """Get the coordinates of an openbabel molecule."""
    coords = np.zeros((obmol.NumAtoms(), 3), float)
    for i, atom in enumerate(openbabel.OBMolAtomIter(obmol)):
        coords[i][0] = atom.x()
        coords[i][1] = atom.y()
        coords[i][2] = atom.z()
    return coords


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
    for r in openbabel.OBMolRingIter(obmol):
        atoms = [x - 1 for x in r._path]  # convert Idx to index
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
        piece = loop[i:i + 4]
        piece2 = piece[::-1]
        fours.append(tuple(piece))
        fours.append(tuple(piece2))
    return fours


def ring_bonds(ring):
    """Get all pairs of atom indices (a,b) where a and b are in the same ring.

    Parameters
    ----------
    ring : tuple(int, bool, str, list(int))
        A ring contains a:
        [0] size (num atoms in ring).
        [1] True = is aromatic, False = not aromatic
        [2] name of ring type (ex: benzene)
        [3] atom indices in connected order around ring

    Returns
    -------
    pairs : list(tuple(int,int))
        All sets of connected atoms in the ring (incl.
        clockwise and counter-clockwise directions).
        Note both indices must be in the ring.

    """
    ring_size = ring[0]
    loop = ring[3] + ring[3][0:1]
    pairs = []
    for i in range(ring_size):
        piece = loop[i:i + 2]
        piece2 = piece[::-1]
        pairs.append(tuple(piece))
        pairs.append(tuple(piece2))
    return pairs


def remove_ring_dihed(rings, diheds):
    """Remove dihedral angles that are only in a ring.

    Parameters
    ----------
    rings : list(int, bool, str, list(int))
        As generated from get_rings function in this module.
    diheds : list(tuple(int,int,int,int))
        Each element is a tuple representing the 4 atom
        indices of the dihedral angle.

    Notes
    -----
    Dihedrals have the form a,b,c,d, where the letters
    represent atom indices. If the b-c bond is within
    a ring, then this dihedral is considered to be a
    ring dihedral, and is removed by this function.
    The c-b case is also checked to account for symmetry.

    Returns
    -------
    diheds : list(tuple(int,int,int,int))
        Updated set of dihedral angles, except with
        the ring dihedral angles (defined in Notes)
        removed.

    """
    to_remove = []
    for ring in rings:
        pairs = ring_bonds(ring)
        for dihed in diheds:
            if dihed[1:3] in pairs:
                to_remove.append(dihed)
    # make sure to only try and remove the dihedral
    # angle once, in case it is part of two rings
    to_remove = set(to_remove)
    for remove in to_remove:
        diheds.remove(remove)
    return diheds


def get_diheds(xyz_file, charge, multip, select_min=True, include_vals=False):
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
    select_min : bool
        Defaults to True, which means select
        a minimum set of dihedral angles (one
        per rotatable bond). False means select
        all dihedral angles.
    include_vals : bool
        Defaults to False, which means the indicies
        of atoms involved is returned. If True, the values
        of the dihedral angles are also returned in the
        tuple.

    Raises
    ------
    AssertionError
        xyz_file does not exist or is not a file.

    Returns
    -------
    If include_vals is False:
    diheds : list(tuple(int,int,int,int))
        Tuple containing a tuple of 4 atom indices. The
        dihedral a-b-c-d as a tuple is (a,b,c,d).
    If include_vals is True:
    diheds : list(tuple(int,int,int,int,float))
        Tuple containing a tuple of 4 atom indices and one
        dihedral angle size in radians. The
        dihedral a-b-c-d as a tuple is (a,b,c,d,size).

    """
    mol = construct_internal(xyz_file, charge, multip, select_min)
    diheds = []
    if not include_vals:
        for intern_coord in mol.ic:
            if isinstance(intern_coord, DihedralAngle):
                atoms = intern_coord.atoms
                # GOpt atoms are np.int64 (not sure openbabel will like that)
                diheds.append((
                    int(atoms[0]), int(atoms[1]),
                    int(atoms[2]), int(atoms[3])
                ))
    else:
        for intern_coord in mol.ic:
            if isinstance(intern_coord, DihedralAngle):
                atoms = intern_coord.atoms
                diheds.append((
                    int(atoms[0]), int(atoms[1]),
                    int(atoms[2]), int(atoms[3]),
                    float(intern_coord.value)
                ))
    return diheds


def construct_internal(xyz_file, charge, multip, select_min=True):
    """Construct an Internal object from an xyzfile."""
    assert os.path.isfile(xyz_file)
    mol = Internal.from_file(xyz_file, charge, multip)
    mol.auto_select_ic(minimum=select_min)
    return mol


def get_new_coordinates(iobj, dihedrals, new_values):
    """Get coordinates from internal object after applying new dihedral angles.

    Parameters
    ----------
    iobj : saddle.internal.Internal object
        An object containing internal coordinate
        information from the original molecular
        coordinates.
    dihedrals : list(tuple(int, int, int, int))
        A list of size-four tuples, where each
        tuple is a set of four indicies (representing
        the connectivity for a dihedral angle).
    new_values : list(len(dihedrals))
        A list of floating point values representing
        the new dihdral values to apply to the
        original geometry. Each index should
        correspond to the tuple in the dihdrals
        list.

    Returns
    -------
    new coordinates as np.array((num_atoms, 3), float)
    or None (if failed to converge)

    """
    # parameters for optimising the geometry
    # increase max_iters or change opt_method if geometry
    # repeatedly fails to converge
    max_iters = 250
    opt_method = "BFGS"
    # make a deepcopy of the saddle object so the original
    # coordinates don't change
    mol = deepcopy(iobj)
    new_ic_values = mol.ic_values
    for i, ic in enumerate(mol.ic):
        try:
            # ic.atoms returns (a,b,c,d) for dihedrals
            # where a-d are integers
            atoms = dihedrals.index(ic.atoms)
        # not a dihedral angle
        except ValueError:
            continue
        # only change it if different
        value = new_values[atoms]
        if ic.value != value:
            new_ic_values[i] = value
    mol.set_target_ic(new_ic_values)
    try:
        mol.optimize_to_target_ic(method=opt_method, max_iter=max_iters)
    except NotConvergeError as e:
        print(f"Warning. Saddle did not converge after {max_iters} iterations of {opt_method}.")
        print(e)
        return None
    return mol.coordinates


def write_coords(coords, atomic_nums, outfile):
    """Write coordinates to an xyz file.

    Parameters
    ----------
    coords : np.array((num_atoms, 3), float)
        The coordinates to write as a numpy array.
    atomic_nums : list(int)
        The atomic numbers that correspond to
        the coordinates, as a list of integers.
    outfile : str
        The name of the output file to write.
        If no directory is given, the current
        working directory is used.

    Notes
    -----
    If the outfile already exists it will be
    overwritten. The comments section of the
    xyz file will be the same as the name of the
    file (minus the extension).

    Returns
    -------
    new_coords : list([char, float, float, float])
        The coordinates that were wrote to the
        outfile.

    """
    file_dir = os.path.dirname(outfile)
    if file_dir == "":
        file_dir = os.getcwd()
        outfile = os.path.join(file_dir, outfile)
    else:
        assert os.path.isdir(file_dir)

    num_atoms = len(coords)
    assert num_atoms == len(atomic_nums)

    comments = os.path.basename(os.path.splitext(os.path.abspath(outfile))[0])
    out_coords = []
    with open(outfile, "w") as f:
        f.write(f"{num_atoms}\n")
        f.write(f"{comments}\n")
        for i, atom in enumerate(coords):
            line = [vetee.tools.periodic_table(atomic_nums[i])] + list(atom)
            out_coords.append(line)
            f.write(f"{' '.join(str(i) for i in line)}\n")
    return out_coords


def get_angles(obmol):
    """Get the angles of an openbabel molecule object.

    Notes
    -----
    Here is the c++ example from the docs for angle iter:
        b = _mol.GetAtom((*angle)[0] + 1);
        a = _mol.GetAtom((*angle)[1] + 1);
        c = _mol.GetAtom((*angle)[2] + 1);
        ang = a->GetAngle(b->GetIdx(), c->GetIdx());
    So each angle from angle iter is of the form (b,a,c).
    The atom has the method, GetAngle, which needs the other
    two connected atoms. Plus 1 for the 1-based
    indexing of atoms.

    The atom.GetAngle(b, c) works where a-b-c and
    b is the vertex. Essentially we get the vectors
    a-b and c-b and find the angle between them.

    Returns
    -------
    a list of length-4 tuples, where each tuple
    has the form: (int, int, int, float)
    The three integers are the indices of
    the atoms that make up the angle (a,b,c, where
    b is the vertex), and the floating point
    value is the angle size in radians.


    """
    angles = []
    for angle in openbabel.OBMolAngleIter(obmol):
        b = obmol.GetAtom(angle[0] + 1)
        a = obmol.GetAtom(angle[1] + 1)
        c = obmol.GetAtom(angle[2] + 1)
        # instead use this method, which takes arguments (a,b,c)
        # returns the angle (in degrees) between the three atoms
        # a, b and c (where a-> b (vertex) -> c )
        value = obmol.GetAngle(a, b, c)
        # convert to radians
        value = geometry_units["degrees"]["radians"](value)
        angles.append((angle[1], angle[0], angle[2], value))
    return angles


def get_torsions(obmol):
    """Get the torsion angles from an Openbabel molecule.

    Parameters
    ----------
    obmol : openbabel.OBMol object
        The molecule for which to collect torsion angles.

    Returns
    -------
    torsions : list(tuples)
        Each tuple is (a,b,c,d,value)
        where a,b,c,d are int indices for atoms
        and value is float of the torsion angle
        in radians.

    """
    torsions = []
    for torsion in openbabel.OBMolTorsionIter(obmol):
        value = get_torsion(obmol, torsion)
        torsions.append((torsion[0], torsion[1], torsion[2], torsion[3], value))
    return torsions
