"""Test the ring module from Kaplan."""

import os

from numpy.testing import assert_raises, assert_allclose

from kaplan.ring import Ring, RingEmptyError, RingOverflowError
from kaplan.pmem import Pmem
from kaplan.inputs import Inputs
from kaplan.tools import TEST_DIR


def test_ring():
    """Test the Ring object from the ring module."""
    inputs = Inputs()
    input_dict = {
        "num_slots": 20,
        "num_geoms": 3,
        "num_atoms": 24,
        "mating_rad": 2,
        "struct_input": os.path.join(TEST_DIR, "caffeine.xyz"),
        "struct_type": "xyz",
        "charge": 0,
        "multip": 1,
        "init_popsize": 15
    }
    inputs.update_inputs(input_dict)
    ring = Ring(20, 15)
    # check that ring has its attributes
    assert inputs.parser.charge == 0
    assert inputs.parser.multip == 1
    assert inputs.num_slots == 20
    assert inputs.init_popsize == 15
    assert ring.num_filled == 15
    assert inputs.num_geoms == 3
    assert len(inputs.coords) == 24
    assert inputs.mating_rad == 2
    assert inputs.fit_form == 0
    assert inputs.coef_energy == 0.5
    assert inputs.coef_rmsd == 0.5
    assert ring.pmems.shape == (20,)
    assert all([ring.pmems[i] is None for i in range(15, 20)])
    # check that the zmatrix is the same as the
    # zmatrix that comes out of the openbabel gui
    czmat = CAFFEINE_ZMATRIX.split('\n')
    ring_czmat = ring.zmatrix.split('\n')
    assert len(czmat) == len(ring_czmat)
    for i, val in enumerate(czmat):
        assert val == ring_czmat[i]


def test_ring_properties():
    inputs = Inputs()
    inputs.update_inputs({
        "struct_input": "propane",
        "init_popsize": 10,
        "num_slots": 50,
    })
    r = Ring(50, 10)
    assert r.occupied == [i for i in range(10)]
    assert_allclose(r.median_energy, -117, atol=2)
    assert_allclose(r.mean_energy, -117, atol=2)
    assert_allclose(r.median_rmsd, 1, atol=0.2)
    assert_allclose(r.mean_rmsd, 1, atol=0.2)
    assert_allclose(r.median_fitness, 297, atol=5)
    assert_allclose(r.mean_fitness, 297, atol=5)
    assert_allclose(r.stdev_fitness, 0.5, atol=0.25)
    assert_allclose(r.stdev_energy, 0.003, atol=0.003)
    assert_allclose(r.stdev_rmsd, 0.26, atol=0.05)



def test_ring_fill():
    """Test the Ring.fill method."""
    inputs = Inputs()
    inputs.update_inputs({
        "struct_input": os.path.join(TEST_DIR, "1,3-butadiene.xyz"),
        "struct_type": "xyz",
        "charge": 0,
        "multip": 1,
        "num_slots": 100,
        "init_popsize": 10,
        "num_geoms": 3,
        "num_muts": 0,
        "num_swaps": 2
    })
    ring = Ring(100, 10)
    # try to put more pmems in than there are slots
    assert_raises(RingOverflowError, ring.fill, 200, 0)
    assert inputs.init_popsize == 10
    ring.fill(1, 0)
    assert ring.num_filled == 11
    # same test twice except getitem syntax can be used
    assert ring.pmems[0].birthday == 0
    assert ring[0].birthday == 0
    ring.fill(5, 1)
    assert ring.num_filled == 16
    assert ring[11].birthday == 1
    # test adding without contiguous segment present
    ring[19] = Pmem(19, 3, 10, 2)
    ring.fill(5, 3)
    assert ring.num_filled == 22
    assert ring[19] is not None
    assert ring[22] is None
    assert ring[18].birthday == 3


def test_ring_getitem():
    """Test the Ring.__getitem__ method."""
    inputs = Inputs()
    inputs.update_inputs({
        "struct_input": os.path.join(TEST_DIR, "1,3-butadiene.xyz"),
        "struct_type": "xyz",
        "charge": 0,
        "multip": 1,
        "num_slots": 100,
        "init_popsize": 1,
        "num_geoms": 3,
        "num_muts": 0,
        "num_swaps": 2
    })
    ring = Ring(1, 100)
    ring.fill(5, 0)
    for i in range(6):
        assert ring[i].birthday == 0
    assert ring[16] is None
    with assert_raises(KeyError):
        ring[1000]
        ring["one"]


def test_ring_setitem():
    """Test the Ring.__setitem__ method."""
    inputs = Inputs()
    inputs.update_inputs({
        "struct_input": os.path.join(TEST_DIR, "1,3-butadiene.xyz"),
        "struct_type": "xyz",
        "charge": 0,
        "multip": 1,
        "num_slots": 100,
        "init_popsize": 1,
        "num_geoms": 3,
        "num_muts": 0,
        "num_swaps": 2
    })
    ring = Ring(1, 10)
    with assert_raises(KeyError):
        ring[1000] = Pmem(1000, 3, 10, 1)
        ring["two"] = Pmem(2, 3, 10, 1)
        ring[1] = "string"
    with assert_raises(AssertionError):
        ring[4] = Pmem(3, 3, 10, 1)
        ring[3] = Pmem(4, 2, 10, 1)
        ring[2] = Pmem(2, 3, 9, 1)
    ring.fill(1, 0)
    assert ring.num_filled == 11
    ring[0] = None
    assert ring.num_filled == 10


def test_ring_update():
    """Test the Ring.update method."""
    # parent_index, child, current_mev
    inputs = Inputs()
    inputs.update_inputs({
        "struct_input": os.path.join(TEST_DIR, "1,3-butadiene.xyz"),
        "struct_type": "xyz",
        "charge": 0,
        "multip": 1,
        "num_slots": 100,
        "num_geoms": 3,
        "num_muts": 0,
        "num_swaps": 2,
        "mating_rad": 2,
        "init_popsize": 1,
    })
    ring = Ring(100, 1)
    # ring is empty, try adding random pmems
    # testing backflow
    # as of 28-05-19 1,3-butadiene has 3 dihedral angles (minimum selection)
    pmem_dihedrals = np.array([[0.25, 0.5, 0.75], [0.2, 0.5, 0.2],
                    [0.8, 0.3, 0.7]])
    ring.update(pmem_dihedrals, 0, 0)
    assert ring[0].dihedrals == pmem_dihedrals
    """
    # possible places the update took place
    slots = [0, 1, 2, 13, 14]
    not_slots = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    assert sum(ring[i] is not None for i in slots) == 1
    assert all(ring[i] is None for i in not_slots)
    ring[0] = None
    ring[1] = None
    ring[2] = None
    ring[13] = None
    ring[14] = None
    # test case where mating_rad is zero
    inputs.mating_rad = 0
    ring.update(0, [[1, 2, 3, 4, 5, 6, 7], [1, 2, 3, 4, 2, 1, 1],
                    [7, 6, 5, 4, 3, 2, 1]], 0)
    assert ring[0] is not None
    assert sum(ring[i] is not None for i in range(inputs.num_slots)) == 1
    # test overflow
    ring[0] = None
    inputs.mating_rad = 4
    ring.update(13, [[1, 2, 3, 4, 5, 6, 7], [1, 2, 3, 4, 2, 1, 1],
                     [7, 6, 5, 4, 3, 2, 1]], 0)
    slots = [13, 14, 0, 1, 2, 9, 10, 11, 12]
    assert (sum(ring[i] is not None for i in slots)) == 1
    not_slots = [3, 4, 5, 6, 7, 8]
    assert all(ring[i] is None for i in not_slots)
    # test no overflow or backflow (now, pmem dist is 4)
    for slot in slots:
        ring[slot] = None
    ring.update(7, [[1, 2, 3, 4, 5, 6, 7], [1, 2, 3, 4, 2, 1, 1],
                    [7, 6, 5, 4, 3, 2, 1]], 0)
    slots = [7, 8, 9, 10, 11, 3, 4, 5, 6]
    assert (sum(ring[i] is not None for i in slots)) == 1
    not_slots = [0, 1, 2, 12, 13, 14]
    assert all(ring[i] is None for i in not_slots)
    for slot in slots:
        ring[slot] = None
    # now check that slot is not updated if child has worse fitness
    inputs.mating_rad = 0
    # fitness = 230.09933808553276
    ring[0] = Pmem(0, 3, 10, 0, [[239, 278, 5, 248, 40, 67, 299],
                                 [36, 123, 295, 111, 322, 267, 170],
                                 [61, 130, 26, 139, 290, 238, 331]])
    ring.set_fitness(0)
    # fitness = 77.4576053229711
    ring.update(0, [[132, 272, 40, 226, 44, 154, 339],
                    [182, 119, 106, 157, 194, 244, 168],
                    [95, 81, 202, 261, 197, 166, 161]], 1)
    assert ring[0].dihedrals == [[239, 278, 5, 248, 40, 67, 299],
                                 [36, 123, 295, 111, 322, 267, 170],
                                 [61, 130, 26, 139, 290, 238, 331]]
    assert ring.num_filled == 1
    assert ring[0].birthday == 0
    """


CAFFEINE_ZMATRIX = """#Put Keywords Here, check Charge and Multiplicity.

 caffeine from pubchem

0  1
O
O  1  r2
N  2  r3  1  a3
N  1  r4  2  a4  3  d4
N  3  r5  1  a5  2  d5
N  4  r6  1  a6  2  d6
C  4  r7  1  a7  2  d7
C  7  r8  4  a8  1  d8
C  1  r9  2  a9  3  d9
C  2  r10  1  a10  3  d10
C  6  r11  4  a11  1  d11
C  3  r12  2  a12  1  d12
C  4  r13  1  a13  2  d13
C  5  r14  2  a14  1  d14
H  11  r15  6  a15  4  d15
H  12  r16  3  a16  2  d16
H  12  r17  3  a17  2  d17
H  12  r18  3  a18  2  d18
H  13  r19  4  a19  1  d19
H  13  r20  4  a20  1  d20
H  13  r21  4  a21  1  d21
H  14  r22  5  a22  2  d22
H  14  r23  5  a23  2  d23
H  14  r24  5  a24  2  d24
Variables:
r2= 4.6919
r3= 2.3268
a3=  61.87
r4= 2.9916
a4=  85.81
d4=   0.02
r5= 2.4221
a5=  29.34
d5=   0.01
r6= 2.2293
a6= 123.04
d6= 359.96
r7= 1.3654
a7=  49.28
d7=   0.03
r8= 1.3682
a8= 105.05
d8= 180.06
r9= 1.2281
a9=  29.10
d9=   0.10
r10= 1.2363
a10=  30.91
d10= 359.97
r11= 1.3171
a11=  34.75
d11= 180.06
r12= 1.4577
a12=  93.57
d12= 179.96
r13= 1.4434
a13=  78.22
d13= 180.01
r14= 1.4593
a14=  89.99
d14= 180.59
r15= 1.0816
a15= 125.45
d15= 180.00
r16= 1.0944
a16= 108.92
d16= 120.20
r17= 1.0929
a17= 111.69
d17= 359.99
r18= 1.0944
a18= 108.91
d18= 239.78
r19= 1.0930
a19= 109.52
d19= 180.00
r20= 1.0927
a20= 108.47
d20= 299.76
r21= 1.0928
a21= 108.46
d21=  60.24
r22= 1.0922
a22= 112.87
d22= 180.08
r23= 1.0947
a23= 108.58
d23=  58.77
r24= 1.0947
a24= 108.59
d24= 301.39

"""

test_ring_properties()