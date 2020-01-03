"""

This module provides various ways to interact
with the Openbabel molecule object, such as
construction from file or string, iteration
over rings, collection of angles or torsions,
and getting and setting of coordinates.

It also contains a very short list of geometry
units for easy conversion. It is independent of the
other Kaplan modules, including inputs.

"""

import openbabel
import pybel
import numpy as np
import os

from copy import deepcopy


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


def create_obmol(obmol_input, input_type, is_file):
    """Create an Openbabel molecule object.

    Parameters
    ----------
    obmol_input : str
        Input value - could be a path+filename or
        a SMILES string, etc.
    input_type : str
        Should be one of:
        com, xyz, inchi, smiles, sdf
        Or another supported file format.
    is_file : bool
        True if obmol_input is a path+filename.
        False if obmol_input is a string.

    Notes
    -----
    Openbabel cannot read com file as input unless
    it is z-matrix format (gzmat is the input_type).
    If an sdf file contains more than one molecule
    or geometry, only the first one will be read.
    Openbabel will not raise an error if the number
    of atoms does not agree with the number of atoms
    listed in an xyz file. Instead, it only reads
    the number of indicated atoms from the input file
    and the resulting object may have a non-standard
    total spin multiplicity.
    This function does not do any sort of geometry
    optimisation, which is recommended if building
    from a smiles or inchi string.

    Returns
    -------
    openbabel.OBMol instance

    """
    assert isinstance(is_file, bool)
    if is_file and not os.path.isfile(obmol_input):
        raise GeometryError(f"No such file: {obmol_input}")

    # create a conversion object instance
    obconvert = openbabel.OBConversion()
    good_format = obconvert.SetInFormat(input_type)
    if not good_format:
        raise GeometryError(f"Unsupported format: {input_type}")
    obmol = openbabel.OBMol()
    if is_file:
        has_mol = obconvert.ReadFile(obmol, obmol_input)
    else:
        has_mol = obconvert.ReadString(obmol, obmol_input)
    if not has_mol:
        raise GeometryError(
            f"Unable to read {input_type} {'file' if is_file else 'string'} input: {obmol_input}")

    requires_builder = ["smiles", "inchi", "smi"]
    if input_type in requires_builder:
        # builder required for identifier-based input
        builder = openbabel.OBBuilder()
        build_success = builder.Build(obmol)
        # unclear how to raise this error, but the builder
        # should return 0 if issues are encountered
        if not build_success:
            raise GeometryError("Could not construct 3D coordinates.")
        obmol.AddHydrogens()

    return obmol


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


def write_coords(coords, atomic_nums, outfile, comments=None):
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
    comments : str
        The comments to put in the xyz file.
        If None (default) writes the name
        of the file to the comments.

    Notes
    -----
    If the outfile already exists it will be
    overwritten. The comments section of the
    xyz file will be the same as the name of the
    file (minus the extension).

    Returns
    -------
    new_coords : list([char, float, float, float])
        The coordinates that were written to the
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

    if comments is None:
        comments = os.path.basename(
            os.path.splitext(os.path.abspath(outfile))[0]
        )
    else:
        comments = comments.replace("\n", " ")
    out_coords = []
    with open(outfile, "w") as f:
        f.write(f"{num_atoms}\n")
        f.write(f"{comments}\n")
        for i, atom in enumerate(coords):
            line = [periodic_table(atomic_nums[i])] + list(atom)
            out_coords.append(line)
            f.write(f"{line[0]:<3} {line[1]:>20f} {line[2]:>20f} {line[3]:>20f}\n")

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


def filter_duplicate_diheds(dihedrals_list, atomic_nums):
    """Remove duplicate dihedrals based on b-c rotatable bond.

    Parameters
    ----------
    dihedrals_list : list(tuple(int, int, int, int, float))
        All possible dihedrals for a molecule.
    atomic_nums : list(int)
        Atomic numbers in order of index, 0-based.
        For example: [1,1,6] means 0 H, 1 H, 2 C.

    Returns
    -------
    list(tuple(int, int, int, int))
        Unique dihedrals with priority given to
        highest atomic number sum.

    """
    atomic_num_sums = []
    for i, dihed in enumerate(dihedrals_list):
        atomic_num_sums.append((sum([atomic_nums[dihed[0]], atomic_nums[dihed[1]],
                                     atomic_nums[dihed[2]], atomic_nums[dihed[3]]]), i))
    atomic_num_sums = sorted(atomic_num_sums, reverse=True)

    dihedrals_by_mass = []
    for i, dihed in enumerate(dihedrals_list):
        dihedrals_by_mass.append(dihedrals_list[atomic_num_sums[i][1]])

    min_diheds_new = []
    unique_bc = []
    for dihed in dihedrals_by_mass:
        if dihed[1:3] not in unique_bc and dihed[2:0:-1] not in unique_bc:
            unique_bc.append(dihed[1:3])
            min_diheds_new.append(dihed[:4])
    return min_diheds_new


class PTableError(Exception):
    """Raised when periodic table cannot generate values."""


def periodic_table(query):
    """Convert atomic numbers to element symbols or vice-versa.

    Parameters
    ----------
    query : list
        If list(int) - atomic numbers.
        If list(str) - element symbols.
    query : int
        atomic number
    query : str
        atomic symbol

    Raises
    ------
    IndexError:
        If given atomic number doesn't exist (needs to
        be in range of 1-118).
    TypeError:
        Input list contains more than one type of input.
        For example, having both str and int in the list.
    ValueError:
        The input str for one of the atoms does not exist
        in the periodic table.

    Returns
    -------
    A list of element symbols (str) (if query was
    type list(int)).
    A list of atomic numbers (int) (if query was
    type list(str)).
    An int (if query was str).
    A str (if query was int).

    """
    atoms = [
        "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
        "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
        "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe",
        "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
        "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
        "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
        "Cs", "Ba",
        "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
        "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",
        "Fr", "Ra",
        "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
        "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
    ]

    try:
        if isinstance(query, str):
            return atoms.index(query) + 1
        elif isinstance(query, int):
            return atoms[query - 1]
        return [atoms[i - 1] for i in query]
    except IndexError:
        raise PTableError("Atomic numbers should be in the range 1-118.") from None
    except TypeError:
        if any([type(a) is not str for a in query]):
            raise PTableError("The input_list should be all integers or all strings.") from None
        try:
            return [atoms.index(i) + 1 for i in query]
        except ValueError:
            raise PTableError("One of the input atoms is not in the periodic table.") from None
    # query is not in list (atoms)
    except ValueError:
        raise PTableError(f"No such element in the periodic table: {query}") from None
