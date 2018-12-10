
import os

from vetee.xyz import Xyz
from numpy.testing import assert_raises

from kaplan.ring import Ring, RingEmptyError, RingOverflowError
from kaplan.pmem import Pmem

# directory for this test file
test_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'testfiles')
 
def test_Ring():
    num_slots = 10
    num_filled = 2
    num_geoms = 3
    num_atoms = 24
    t_size = 7
    num_muts = 3
    num_swaps = 1
    pmem_dist = 2
    parser = Xyz(os.path.join(test_dir, "caffeine.xyz"))
    parser.charge = 0
    parser.multip = 1
    ring = Ring(num_geoms, num_atoms, num_slots, pmem_dist, 0, 0.5, 0.5, parser)
    # check that ring has its attributes
    assert ring.parser.charge == 0
    assert ring.parser.multip == 1
    assert ring.num_slots == 10
    assert ring.num_filled == 0
    assert ring.num_geoms == 3
    assert ring.num_atoms == 24
    assert ring.pmem_dist == 2
    assert ring.fit_form == 0
    assert ring.coef_energy == 0.5
    assert ring.coef_rmsd == 0.5
    assert ring.pmems.shape == (10,)
    assert all([ring.pmems[i] is None for i in range(10)])
    # check that the zmatrix is the same as the
    # zmatrix that comes out of the openbabel gui
    czmat = caffeine_zmatrix.split('\n')
    ring_czmat = ring.zmatrix.split('\n')
    assert len(czmat) == len(ring_czmat)
    for i in range(len(czmat)):
        assert czmat[i] == ring_czmat[i]

def test_Ring_fill():
    parser = Xyz(os.path.join(test_dir, "1,3-butadiene.xyz"))
    parser.charge = 0
    parser.multip = 1
    ring = Ring(3, 10, 15, 2, 0, 0.5, 0.5, parser)
    # try to put more pmems in than there are slots
    assert_raises(RingOverflowError, ring.fill, 200, 0)
    assert ring.num_filled == 0
    ring.fill(1, 0)
    assert ring.num_filled == 1
    # same test twice except getitem syntax can be used
    assert ring.pmems[0].birthday == 0
    assert ring[0].birthday == 0
    ring.fill(5, 1)
    assert ring.num_filled == 6
    assert ring[1].birthday == 1
    # test adding without contiguous segment present
    ring[9] = Pmem(9, 3, 10, 2)
    ring.fill(5, 3)
    assert ring.num_filled == 12
    assert ring[11] is not None
    assert ring[12] is None
    assert ring[8].birthday == 3

def test_Ring_getitem():
    parser = Xyz(os.path.join(test_dir, "1,3-butadiene.xyz"))
    parser.charge = 0
    parser.multip = 1
    ring = Ring(3, 10, 15, 2, 0, 0.5, 0.5, parser)
    ring.fill(5, 0)
    for i in range(5):
        assert ring[i].birthday == 0
    assert ring[5] is None
    with assert_raises(KeyError):
        ring[1000]
        ring['one']

def test_Ring_setitem():
    parser = Xyz(os.path.join(test_dir, "1,3-butadiene.xyz"))
    parser.charge = 0
    parser.multip = 1
    ring = Ring(3, 10, 15, 2, 0, 0.5, 0.5, parser)
    with assert_raises(KeyError):
        ring[1000] = Pmem(1000, 3, 10, 1)
        ring['two'] = Pmem(2, 3, 10, 1)
        ring[1] = 'string'
    with assert_raises(AssertionError):
        ring[4] = Pmem(3, 3, 10, 1)
        ring[3] = Pmem(4, 2, 10, 1)
        ring[2] = Pmem(2, 3, 9, 1)
    ring.fill(1,0)
    assert ring.num_filled == 1
    ring[0] = None
    assert ring.num_filled == 0

def test_Ring_update():
    # parent_index, child, current_mev
    parser = Xyz(os.path.join(test_dir, "1,3-butadiene.xyz"))
    parser.charge = 0
    parser.multip = 1
    ring = Ring(3, 10, 15, 2, 0, 0.5, 0.5, parser)
    # ring is empty, try adding random pmems
    # testing backflow
    ring.update(0, [[1,2,3,4,5,6,7],[1,2,3,4,2,1,1],[7,6,5,4,3,2,1]], 0)
    # possible places the update took place
    slots = [0,1,2,13,14]
    not_slots = [3,4,5,6,7,8,9,10,11,12]
    assert sum(ring[i] is not None for i in slots) == 1
    assert all(ring[i] is None for i in not_slots)
    ring[0] = None
    ring[1] = None
    ring[2] = None
    ring[13] = None
    ring[14] = None
    # test case where pmem_dist is zero
    ring.pmem_dist = 0
    ring.update(0, [[1,2,3,4,5,6,7],[1,2,3,4,2,1,1],[7,6,5,4,3,2,1]], 0)
    assert ring[0] is not None
    assert sum(ring[i] is not None for i in range(ring.num_slots)) == 1
    # test overflow
    ring[0] = None
    ring.pmem_dist = 4
    ring.update(13, [[1,2,3,4,5,6,7],[1,2,3,4,2,1,1],[7,6,5,4,3,2,1]], 0)
    slots = [13, 14, 0, 1, 2, 9, 10, 11, 12]
    assert (sum(ring[i] is not None for i in slots)) == 1
    not_slots = [3,4,5,6,7,8]
    assert all(ring[i] is None for i in not_slots)
    # test no overflow or backflow (now, pmem dist is 4)
    for slot in slots:
        ring[slot] = None
    ring.update(7, [[1,2,3,4,5,6,7],[1,2,3,4,2,1,1],[7,6,5,4,3,2,1]], 0)
    slots = [7, 8, 9, 10, 11, 3, 4, 5, 6]
    assert (sum(ring[i] is not None for i in slots)) == 1
    not_slots = [0,1,2,12,13,14]
    assert all(ring[i] is None for i in not_slots)
    for slot in slots:
        ring[slot] = None
    # now check that slot is not updated if child has worse fitness
    ring.pmem_dist = 0
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

caffeine_zmatrix = """#Put Keywords Here, check Charge and Multiplicity.

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

if __name__ == "__main__":
    test_Ring()
    test_Ring_fill()
    test_Ring_getitem()
