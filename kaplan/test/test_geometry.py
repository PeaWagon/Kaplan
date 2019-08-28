"""Test geometry module from Kaplan."""

import os
import openbabel
import pybel
import numpy as np

from math import pi

from numpy.testing import assert_raises
from vetee.coordinates import CoordinatesError
from vetee.tools import periodic_table

from kaplan.geometry import get_diheds, create_obmol, get_rings,\
    update_obmol, remove_ring_dihed, construct_fours,\
    get_struct_info, geometry_units, get_coords
from kaplan.inputs import Inputs
from kaplan.tools import TEST_DIR, amino_acids


def test_get_diheds():
    test1 = os.path.join(TEST_DIR, "2-pentanol.xyz")
    diheds = get_diheds(test1, 0, 1)
    assert len(diheds) == 5
    for i in range(5):
        assert any([
            (17, 0, 2, 4) == diheds[i],
            (3, 1, 2, 4) == diheds[i],
            (2, 1, 3, 5) == diheds[i],
            (1, 2, 4, 13) == diheds[i],
            (1, 3, 5, 16) == diheds[i],
        ])
        assert isinstance(diheds[i], tuple)
        assert isinstance(diheds[i][0], int)
        assert isinstance(diheds[i][1], int)
        assert isinstance(diheds[i][2], int)
        assert isinstance(diheds[i][3], int)


def test_create_obmol():
    test1 = os.path.join(TEST_DIR, "2-pentanol.xyz")
    obmol = create_obmol(test1, 0, 1)
    atoms = [
        "O", "C", "C", "C", "C", "C", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H"
    ]
    atom_nums = [8, 6, 6, 6, 6, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    for i, atom in enumerate(openbabel.OBMolAtomIter(obmol)):
        letter_num = atom.GetAtomicNum()
        assert letter_num == atom_nums[i]
        letter = periodic_table(letter_num)
        assert letter == atoms[i]
        assert atom.x() != 0
        assert atom.y() != 0
        assert atom.z() != 0


def test_get_rings():
    test1 = os.path.join(TEST_DIR, "caffeine.xyz")
    obmol = create_obmol(test1, 0, 1)
    # to construct the drawing of caffeine:
    # for bond in openbabel.OBMolBondIter(obmol):
    #    print(bond.GetBeginAtom().GetIdx(), bond.GetEndAtom().GetIdx(), bond.GetBondOrder())
    rings = get_rings(obmol)
    assert rings == [
        (5, True, "imidazole", [5, 7, 6, 3, 10]),
        (6, True, "pyrimidine", [4, 8, 6, 7, 2, 9]),
    ]


def test_construct_fours():
    test1 = os.path.join(TEST_DIR, "caffeine.xyz")
    obmol = create_obmol(test1, 0, 1)
    rings = get_rings(obmol)
    ring1, ring2 = rings  # pylint: disable=unbalanced-tuple-unpacking
    assert len(construct_fours(ring1)) == 10
    assert len(construct_fours(ring2)) == 12
    three_membered_ring = (3, False, "oxirane", [0, 1, 2])
    assert_raises(AssertionError, construct_fours, three_membered_ring)
    four_membered_ring = (4, False, "oxetane", [0, 1, 2, 3])
    result = construct_fours(four_membered_ring)
    assert len(result) == 8
    # clockwise
    assert (0, 1, 2, 3) in result
    assert (1, 2, 3, 0) in result
    assert (2, 3, 0, 1) in result
    assert (3, 0, 1, 2) in result
    # counter-clockwise
    assert (3, 2, 1, 0) in result
    assert (0, 3, 2, 1) in result
    assert (1, 0, 3, 2) in result
    assert (2, 1, 0, 3) in result
    five_membered_ring = (5, False, "tetrahydrofuran", [6, 2, 7, 1, 3])
    six_membered_ring = (6, True, "benzene", [3, 5, 2, 6, 8, 1])
    seven_membered_ring = (7, False, "oxepane", [11, 22, 33, 22, 5, 1, 3])
    result = construct_fours(five_membered_ring)
    assert len(result) == 10
    result = construct_fours(six_membered_ring)
    assert len(result) == 12
    result = construct_fours(seven_membered_ring)
    assert len(result) == 14


def test_update_obmol():
    inputs = Inputs()
    test1 = os.path.join(TEST_DIR, "2-pentanol.xyz")
    obmol = create_obmol(test1, 0, 1)
    dihedrals = get_diheds(test1, 0, 1)
    rings = get_rings(obmol)
    assert rings == []
    dihedrals = remove_ring_dihed(rings, dihedrals)
    assert dihedrals != []
    assert len(dihedrals) == 5
    for i in range(5):
        # create new sets of dihedrals
        new_dihed = np.random.choice(
            inputs.avail_diheds, size=len(dihedrals)
        )
        new_coords = update_obmol(obmol, dihedrals, new_dihed)
        for i, coord in enumerate(new_coords):
            print(obmol.GetAtom(i + 1).GetAtomicNum(), coord[0], coord[1], coord[2])
        print("\n\n")

    test2 = os.path.join(TEST_DIR, "caffeine.xyz")
    obmol = create_obmol(test2, 0, 1)
    dihedrals = get_diheds(test2, 0, 1)
    rings = get_rings(obmol)
    dihedrals = remove_ring_dihed(rings, dihedrals)
    assert (10, 5, 7, 6) not in dihedrals
    assert (7, 5, 10, 3) not in dihedrals
    assert (8, 6, 7, 2) not in dihedrals
    assert (7, 6, 8, 4) not in dihedrals
    # print('new:')
    # for dihed in dihedrals:
    #    print([x+1 for x in dihed])
    """
    test2 = os.path.join(TEST_DIR, "oxazole.xyz")
    obmol = create_obmol(test2)
    dihedrals = get_diheds(test2, 0, 1)
    print(dihedrals)
    rings = get_rings(obmol)
    dihedrals = remove_ring_dihed(rings, dihedrals)
    print(dihedrals)
    """

    """Test to see if update_obmol actually changes the dihedral angle.

    Results
    -------
    The update_obmol does update the Openbabel molecule coordinates
    and torsion values. However, it was discovered that torsion
    angles are absolute (not relative to the rest of the geometry,
    as was intially thought). This result means that it is not
    necessary to reset the obmol coordinates each time to the
    starting geometry - this step is nice because it makes
    sure the geometries are centred about the initial structure
    (and the floating point error is reduced, since each application
    of a torsion angle would incur some error).

    """
    infile = os.path.join(TEST_DIR, "butane.xyz")

    inputs = Inputs()
    inputs.update_inputs({
        "struct_input": infile,
        "struct_type": "xyz",
        "charge": 0,
        "multip": 1,
        "prog": "openbabel",
        "output_dir": TEST_DIR,
    })

    new_angle = 0.0
    torsion_original = inputs.obmol.GetTorsion(3, 1, 2, 4)
    t_orig2 = inputs.obmol.GetTorsion(2, 1, 3, 11)
    t_orig3 = inputs.obmol.GetTorsion(1, 2, 4, 14)
    # print(torsion_original, t_orig2, t_orig3)
    assert np.allclose(
        [torsion_original, t_orig2, t_orig3],
        [geometry_units["radians"]["degrees"](3.1215158160645746),
         geometry_units["radians"]["degrees"](3.1405551846521624),
         geometry_units["radians"]["degrees"](1.0507839153871283)]
    )

    d = inputs.obmol.GetAtom(4)
    result = d.GetVector()
    assert np.allclose([result.GetX(), result.GetY(), result.GetZ()], inputs.coords[3])
    # print(result.GetX(), result.GetY(), result.GetZ())

    torsions = [new_angle, 3.1405551846521624, 1.0507839153871283]
    new_coords = update_obmol(inputs.obmol, inputs.diheds, torsions)
    obmol_coords = get_coords(inputs.obmol)
    # GetTorsion returns degrees, SetTorsion takes radians
    # most Openbabel functions work with IdX, which is index + 1
    updated_torsion = inputs.obmol.GetTorsion(3, 1, 2, 4)
    t_new2 = inputs.obmol.GetTorsion(2, 1, 3, 11)
    t_new3 = inputs.obmol.GetTorsion(1, 2, 4, 14)

    # (2, 0, 1, 3), 178.8496819438351 3.1215158160645746
    # (1, 0, 2, 10), 179.94055740849785 3.1405551846521624
    # (0, 1, 3, 13), 60.205483531914254 1.0507839153871283

    # print(updated_torsion, t_new2, t_new3)
    assert np.allclose(new_coords, obmol_coords)
    assert np.allclose(
        [updated_torsion, t_new2, t_new3],
        geometry_units["radians"]["degrees"](np.array(torsions))
    )

    d = inputs.obmol.GetAtom(4)
    result = d.GetVector()
    assert not np.allclose([result.GetX(), result.GetY(), result.GetZ()], inputs.coords[3])
    # print(result.GetX(), result.GetY(), result.GetZ())
    # print("\n")

    new_angle = pi / 2
    updated_torsion = inputs.obmol.GetTorsion(3, 1, 2, 4)
    t_new2 = inputs.obmol.GetTorsion(2, 1, 3, 11)
    t_new3 = inputs.obmol.GetTorsion(1, 2, 4, 14)
    # print(updated_torsion, t_new2, t_new3)
    torsions = [new_angle, 3.1405551846521624, 1.0507839153871283]
    new_coords = update_obmol(inputs.obmol, inputs.diheds, torsions)
    obmol_coords = get_coords(inputs.obmol)
    updated_torsion = inputs.obmol.GetTorsion(3, 1, 2, 4)
    t_new2 = inputs.obmol.GetTorsion(2, 1, 3, 11)
    t_new3 = inputs.obmol.GetTorsion(1, 2, 4, 14)
    assert np.allclose(
        [updated_torsion, t_new2, t_new3],
        geometry_units["radians"]["degrees"](np.array(torsions))
    )
    # print(updated_torsion, t_new2, t_new3)
    assert np.allclose(new_coords, obmol_coords)

    d = inputs.obmol.GetAtom(4)
    result = d.GetVector()
    assert not np.allclose([result.GetX(), result.GetY(), result.GetZ()], inputs.coords[3])
    # print(result.GetX(), result.GetY(), result.GetZ())
    # print("\n")

    # print(new_coords, "\n")
    # print(obmol_coords)
    torsion_original = inputs.obmol.GetTorsion(3, 1, 2, 4)
    t_orig2 = inputs.obmol.GetTorsion(2, 1, 3, 11)
    t_orig3 = inputs.obmol.GetTorsion(1, 2, 4, 14)
    # print(torsion_original, t_orig2, t_orig3)
    torsions = [0.0, pi, pi / 3]
    new_coords = update_obmol(inputs.obmol, inputs.diheds, torsions)
    obmol_coords = get_coords(inputs.obmol)
    updated_torsion = inputs.obmol.GetTorsion(3, 1, 2, 4)
    t_new2 = inputs.obmol.GetTorsion(2, 1, 3, 11)
    t_new3 = inputs.obmol.GetTorsion(1, 2, 4, 14)
    assert np.allclose(
        [updated_torsion, t_new2, t_new3],
        geometry_units["radians"]["degrees"](np.array(torsions))
    )
    # print(updated_torsion, t_new2, t_new3)
    assert np.allclose(new_coords, obmol_coords)


def test_atom_indices():
    test1 = os.path.join(TEST_DIR, "caffeine.xyz")
    obmol = create_obmol(test1, 0, 1)
    for atom in openbabel.OBMolAtomIter(obmol):
        print(atom.GetIdx(), atom.GetAtomicNum(), atom.x(), atom.y(), atom.z())


def test_get_struct_info():
    test_names = [
        "propane",
        "cysteine",
        "alanine",
        "tryptophan",
        "tyrosine",
        "glycine",
        "threonine",
        "arginine",
        "proline",
        "glutamine",
        "isoleucine",
        "leucine",
    ]
    for test in test_names:
        result = get_struct_info(test)
        assert result.value == test
        assert result.identifier == "name"
        assert result.charge == 0
        assert result.multip == 1
        assert len(result.coords)
    failures = [
        100403,
        "not-a-molecule",
        "CCCC",
        "InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3",  # inchi
    ]
    # note inchikey is accepted as a name, so it won't raise an
    # error (unless inchikey is not in synonyms)
    for f in failures:
        assert_raises(CoordinatesError, get_struct_info, f)
    result = get_struct_info("ammonium")
    assert result.charge == 1
    assert len(result.coords) == 5
    assert result.multip == 1
    assert result.num_atoms == 5
    assert result.atomic_nums == [7, 1, 1, 1, 1]


def test_get_coords():
    # see if obatom.x() and obatom.GetX() return the same values
    for aa in amino_acids:
        sdf_file = os.path.join(TEST_DIR, f"{aa}.sdf")
        assert os.path.isfile(sdf_file)
        mol = pybel.readfile("sdf", sdf_file)
        mol = mol.__next__()
        obmol = mol.OBMol
        coords_from_get = np.zeros((obmol.NumAtoms(), 3), float)
        for i, atom in enumerate(openbabel.OBMolAtomIter(obmol)):
            coords_from_get[i][0] = atom.GetX()
            coords_from_get[i][1] = atom.GetY()
            coords_from_get[i][2] = atom.GetZ()

        # read in the file again just to be sure
        mol2 = pybel.readfile("sdf", sdf_file)
        mol2 = mol2.__next__()
        obmol2 = mol2.OBMol
        coords_from_func = get_coords(obmol2)

        assert np.allclose(coords_from_get, coords_from_func)

    # make sure the coordinates are the same when reading
    # an sdf file
    glycine_file_coords = np.array([
        [-1.6487, 0.6571, -0.0104],
        [-0.4837, -1.2934, -0.0005],
        [1.9006, -0.0812, -0.0090],
        [0.7341, 0.7867, 0.0079],
        [-0.5023, -0.0691, 0.0120],
        [0.7326, 1.4215, -0.8824],
        [0.7464, 1.4088, 0.9069],
        [1.8743, -0.6844, -0.8301],
        [1.8887, -0.6969, 0.8031],
        [-2.4447, 0.0839, -0.0260]
    ])
    mol3 = pybel.readfile("sdf", os.path.join(TEST_DIR, "glycine.sdf"))
    mol3 = mol3.__next__()
    obmol3 = mol3.OBMol
    glycine_mol_coords = get_coords(obmol3)
    assert np.allclose(glycine_file_coords, glycine_mol_coords)
