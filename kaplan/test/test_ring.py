"""Test the ring module from Kaplan."""

import os
import numpy as np

from numpy.testing import assert_raises, assert_allclose

from kaplan.ring import Ring, RingEmptyError, RingOverflowError
from kaplan.pmem import Pmem
from kaplan.inputs import Inputs
from kaplan.tools import TEST_DIR
from kaplan.geometry import geometry_units


def test_ring():
    """Test the Ring object from the ring module."""
    inputs = Inputs()
    input_dict = {
        "num_slots": 20,
        "num_geoms": 3,
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
    assert inputs.charge == 0
    assert inputs.multip == 1
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


def test_ring_properties():
    inputs = Inputs()
    inputs.update_inputs({
        "struct_input": "propane",
        "init_popsize": 10,
        "num_slots": 50,
        "normalise": False,
    })
    r = Ring(50, 10)
    assert r.occupied == [i for i in range(10)]
    assert_allclose(r.median_energy, -117, atol=2)
    assert_allclose(r.mean_energy, -117, atol=2)
    assert_allclose(r.median_rmsd, 1, atol=0.2)
    assert_allclose(r.mean_rmsd, 1, atol=0.2)
    assert_allclose(r.median_fitness, 297, atol=5)
    assert_allclose(r.mean_fitness, 297, atol=5)
    assert_allclose(r.stdev_fitness, 0.5, atol=0.5)
    assert_allclose(r.stdev_energy, 0.003, atol=0.003)
    assert_allclose(r.stdev_rmsd, 0.26, atol=0.1)


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
    ring[19] = Pmem(19, 3, 3, 3)
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
    ring = Ring(100, 1)
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
    ring = Ring(100, 1)
    # slot out of range
    with assert_raises(KeyError):
        ring[1000] = Pmem(1000, 3, 10, 1)
    # slot should be an integer
    with assert_raises(KeyError):
        ring["two"] = Pmem(2, 3, 10, 1)
    # slot value should be None or Pmem
    with assert_raises(KeyError):
        ring[1] = "string"
    # pmem doesn't go to the correct slot
    with assert_raises(AssertionError):
        ring[4] = Pmem(3, 3, 10, 1)
    with assert_raises(AssertionError):
        ring[3] = Pmem(4, 2, 10, 1)
    ring.fill(1, 0)
    assert ring.num_filled == 2
    ring[0] = None
    assert ring.num_filled == 1


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
    pmem_dihedrals = np.array([
        [0.25, 0.5, 0.75], [0.2, 0.5, 0.2], [0.8, 0.3, 0.7]
    ])
    ring.update(pmem_dihedrals, 1, 0)
    assert (ring[1].dihedrals == pmem_dihedrals).all()
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


def test_ring_iteration():
    """Test the ring iterator and __str__ for pmem and ring."""
    pmem = Pmem(0, 0, 3, 2)
    print(pmem)

    inputs = Inputs()
    inputs.update_inputs({
        "struct_type": "name",
        "struct_input": "propane",
        "charge": 0,
        "multip": 1,
        "init_popsize": 5,
        "num_slots": 20
    })
    ring = Ring(20, 5)

    for i, r in enumerate(ring):
        print(i, r)

    print(ring)

    pmem5 = Pmem(5, 1, inputs.num_geoms, inputs.num_diheds)
    ring[5] = pmem5
    for dihedrals in ring[5]:
        # print(dihedrals)
        print([geometry_units["radians"]["degrees"](g) for g in dihedrals])
    print(ring[5])
    print(pmem5)
