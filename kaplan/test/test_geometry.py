"""Test geometry module from Kaplan."""

import os
import openbabel
import numpy as np

from numpy.testing import assert_raises
from vetee.coordinates import CoordinatesError

from kaplan.geometry import get_min_dihed, create_obmol, get_rings,\
    update_obmol, remove_ring_dihed, construct_fours, MIN_VALUE, MAX_VALUE,\
    get_struct_info
from kaplan.inputs import Inputs

# directory for this test file
TEST_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "testfiles")


def test_get_min_dihed():
    test1 = os.path.join(TEST_DIR, "2-pentanol.xyz")
    min_dihed = get_min_dihed(test1, 0, 1)
    assert len(min_dihed) == 5
    for i in range(5):
        assert any([
            (17,0,2,4) == min_dihed[i][0],
            (3,1,2,4) == min_dihed[i][0],
            (2,1,3,5) == min_dihed[i][0],
            (1,2,4,13) == min_dihed[i][0],
            (1,3,5,16) == min_dihed[i][0],
        ])
        assert isinstance(min_dihed[i][0], tuple)
        assert isinstance(min_dihed[i][1], float)


def test_create_obmol():
    test1 = os.path.join(TEST_DIR, "2-pentanol.xyz")
    obmol = create_obmol(test1)
    for atom in openbabel.OBMolAtomIter(obmol):
        #print(atom.x(), atom.y(), atom.z())
        #print(atom.GetAtomicNum())
        pass


def test_get_rings():
    test1 = os.path.join(TEST_DIR, "caffeine.xyz")
    obmol = create_obmol(test1)
    # to construct the drawing of caffeine:
    #for bond in openbabel.OBMolBondIter(obmol):
    #    print(bond.GetBeginAtom().GetIdx(), bond.GetEndAtom().GetIdx(), bond.GetBondOrder())
    rings = get_rings(obmol)
    assert rings == [
        (5, True, "imidazole", [5, 7, 6, 3, 10]),
        (6, True, "pyrimidine", [4, 8, 6, 7, 2, 9]),
    ]


def test_construct_fours():
    test1 = os.path.join(TEST_DIR, "caffeine.xyz")
    obmol = create_obmol(test1)
    rings = get_rings(obmol)
    ring1, ring2 = rings
    assert len(construct_fours(ring1)) == 10
    assert len(construct_fours(ring2)) == 12
    three_membered_ring = (3, False, "oxirane", [0,1,2])
    assert_raises(AssertionError, construct_fours, three_membered_ring)
    four_membered_ring = (4, False, "oxetane", [0,1,2,3])
    result = construct_fours(four_membered_ring)
    assert len(result) == 8
    # clockwise
    assert (0,1,2,3) in result
    assert (1,2,3,0) in result
    assert (2,3,0,1) in result
    assert (3,0,1,2) in result
    # counter-clockwise
    assert (3,2,1,0) in result
    assert (0,3,2,1) in result
    assert (1,0,3,2) in result
    assert (2,1,0,3) in result
    five_membered_ring = (5, False, "tetrahydrofuran", [6,2,7,1,3])
    six_membered_ring = (6, True, "benzene", [3,5,2,6,8,1])
    seven_membered_ring = (7, False, "oxepane", [11,22,33,22,5,1,3])
    result = construct_fours(five_membered_ring)
    assert len(result) == 10
    result = construct_fours(six_membered_ring)
    assert len(result) == 12
    result = construct_fours(seven_membered_ring)
    assert len(result) == 14



def test_update_obmol():
    test1 = os.path.join(TEST_DIR, "2-pentanol.xyz")
    obmol = create_obmol(test1)
    dihedrals = get_min_dihed(test1, 0, 1)
    dihedrals = [x[0] for x in dihedrals] # remove dihedral values
    rings = get_rings(obmol)
    assert rings == []
    dihedrals = remove_ring_dihed(rings, dihedrals)
    assert dihedrals != []
    assert len(dihedrals) == 5
    for i in range(5):
        # create new sets of dihedrals
        new_dihed = np.random.uniform(MIN_VALUE, MAX_VALUE,
                                      size=len(dihedrals))
        new_coords = update_obmol(obmol, dihedrals, new_dihed)
        for i, coord in enumerate(new_coords):
            print(obmol.GetAtom(i+1).GetAtomicNum(), coord[0], coord[1], coord[2])
        print("\n\n")




    test2 = os.path.join(TEST_DIR, "caffeine.xyz")
    obmol = create_obmol(test2)
    dihedrals = get_min_dihed(test2, 0, 1)
    dihedrals = [x[0] for x in dihedrals] # remove dihedral values
    rings = get_rings(obmol)
    dihedrals = remove_ring_dihed(rings, dihedrals)
    assert (10, 5, 7, 6) not in dihedrals
    assert (7, 5, 10, 3) not in dihedrals
    assert (8, 6, 7, 2) not in dihedrals
    assert (7, 6, 8, 4) not in dihedrals
    #print('new:')
    #for dihed in dihedrals:
    #    print([x+1 for x in dihed])
    """
    test2 = os.path.join(TEST_DIR, "oxazole.xyz")
    obmol = create_obmol(test2)
    dihedrals = get_min_dihed(test2, 0, 1)
    dihedrals = [x[0] for x in dihedrals] # remove dihedral values
    print(dihedrals)
    rings = get_rings(obmol)
    dihedrals = remove_ring_dihed(rings, dihedrals)
    print(dihedrals)
    """

def test_atom_indices():
    test1 = os.path.join(TEST_DIR, "caffeine.xyz")
    obmol = create_obmol(test1)
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
        "InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3", # inchi
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
    assert result.atomic_nums == [7,1,1,1,1]
    


#test_get_min_dihed()
#test_create_obmol()
#test_get_rings()
#test_construct_fours()
#test_update_obmol()
#test_atom_indices()
test_get_struct_info()
