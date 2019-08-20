from kaplan.ring import Ring, RingEmptyError
from kaplan.pmem import Pmem
from kaplan.extinction import asteroid, deluge, agathic, plague, apply_extinction
from kaplan.inputs import Inputs, InputError
from random import randint
from statistics import mean
from numpy.testing import assert_raises


def test_apply_extinction():
    inputs = Inputs()
    inputs.num_slots = 100
    inputs.num_geoms = 0
    inputs.num_diheds = 0

    # check that an error is raised if no pmems have fitness
    # values set
    ring1 = Ring(inputs.num_slots, 0)
    for i in range(inputs.num_slots):
        pmem = Pmem(i, 0, 0, inputs.num_diheds)
        ring1[i] = pmem
    assert_raises(RingEmptyError, apply_extinction, ring1, "deluge", inputs.normalise)

    # check that an error is raised for a ring with no pmems
    ring2 = Ring(inputs.num_slots, 0)
    assert_raises(RingEmptyError, apply_extinction, ring2, "deluge", inputs.normalise)

    inputs = Inputs()
    inputs.num_slots = 10

    # if slots are not in range of ring, raise an error
    ring3 = Ring(inputs.num_slots, 0)
    for i in range(inputs.num_slots):
        pmem = Pmem(i, 0, 0, inputs.num_diheds)
        ring3[i] = pmem
        ring3[i]._ring_loc = -4 * (i + inputs.num_slots)
    assert_raises(InputError, apply_extinction, ring3, "asteroid", inputs.normalise)

    ring4 = Ring(inputs.num_slots, 0)
    for i in range(inputs.num_slots):
        pmem = Pmem(i, 0, 0, inputs.num_diheds)
        ring4[i] = pmem
        ring4[i]._ring_loc = i + inputs.num_slots
    assert_raises(InputError, apply_extinction, ring4, "asteroid", inputs.normalise)


def test_plague():
    inputs = Inputs()
    inputs.num_slots = 10
    inputs.num_geoms = 0
    inputs.normalise = False
    inputs.num_diheds = 3

    ring = Ring(inputs.num_slots, 0)
    for i in range(inputs.num_slots):
        pmem = Pmem(i, 0, 0, inputs.num_diheds)
        ring[i] = pmem
        ring[i].fitness = randint(1, 100)
    for i in range(10):
        plague(ring, inputs.normalise)
    assert ring.num_filled >= 1

    # test what happens when all pmems get the same fitness
    inputs.num_slots = 5
    ring2 = Ring(inputs.num_slots, 0)
    for i in range(inputs.num_slots):
        pmem = Pmem(i, 0, 0, inputs.num_diheds)
        ring2[i] = pmem
        ring2[i].fitness = 10
    for i in range(10):
        plague(ring2, inputs.normalise)
    assert ring2.num_filled >= 1


def test_agathic():
    inputs = Inputs()
    inputs.num_slots = 10
    inputs.num_geoms = 0
    inputs.normalise = False
    inputs.num_diheds = 3
    ring = Ring(inputs.num_slots, 0)
    for i in range(inputs.num_slots):
        pmem = Pmem(i, i, 0, inputs.num_diheds)
        ring[i] = pmem
    for i in range(10):
        agathic(ring)
    assert ring.num_filled >= 1

    # test what happens when all pmems get the same bday
    inputs.num_slots = 5
    ring2 = Ring(inputs.num_slots, 0)
    for i in range(inputs.num_slots):
        pmem = Pmem(i, 0, 0, inputs.num_diheds)
        ring2[i] = pmem
    for i in range(10):
        plague(ring2, inputs.normalise)
    assert ring2.num_filled >= 1

    # test with some bdays being None
    ring3 = Ring(inputs.num_slots, 0)
    for i in range(inputs.num_slots):
        pmem = Pmem(i, i, 0, inputs.num_diheds)
        ring3[i] = pmem
    ring3[0]._birthday = None
    ring3[1]._birthday = None
    for i in range(10):
        plague(ring3, inputs.normalise)
    assert ring3.num_filled >= 1


def test_deluge():
    inputs = Inputs()
    inputs.num_slots = 100
    inputs.num_geoms = 0
    inputs.normalise = False
    inputs.num_diheds = 3

    # test that deluge works on random fitness
    # values between 1 and 100
    # make sure that not all pmems got killed
    ring = Ring(inputs.num_slots, 0)
    for i in range(inputs.num_slots):
        pmem = Pmem(i, 0, 0, inputs.num_diheds)
        ring[i] = pmem
        ring[i].fitness = randint(1, 100)
    for i in range(1000):
        deluge(ring, ring.occupied, inputs.normalise)
    assert ring.num_filled > 0
    assert ring.occupied != []
    # very unlikely that any pmems with fitness lower than
    # 90 survivied
    assert mean([ring[i].fitness for i in ring.occupied]) > 90

    # check that the ring can have a deluge
    # even when some fitness values are not initialised
    ring2 = Ring(inputs.num_slots, 0)
    for i in range(inputs.num_slots):
        pmem = Pmem(i, 0, 0, inputs.num_diheds)
        ring2[i] = pmem
        if i < inputs.num_slots / 2:
            ring2[i].fitness = 95
    ring2[-1].fitness = 100
    assert ring2[-1].ring_loc == inputs.num_slots - 1
    for i in range(1000):
        deluge(ring2, ring2.occupied, inputs.normalise)
    assert ring2.num_filled == int(inputs.num_slots / 2) + 1

    # check that there is still a pmem left if
    # the max fitness is zero
    ring5 = Ring(inputs.num_slots, 0)
    for i in range(inputs.num_slots):
        pmem = Pmem(i, 0, 0, inputs.num_diheds)
        ring5[i] = pmem
        ring5[i].fitness = -i
    for i in range(1000):
        deluge(ring5, ring5.occupied, inputs.normalise)
    assert ring5.num_filled == 1
    assert ring5[0].fitness == 0

    # check that a pmem remains when the max
    # fitness is negative
    ring6 = Ring(inputs.num_slots, 0)
    for i in range(inputs.num_slots):
        pmem = Pmem(i, 0, 0, inputs.num_diheds)
        ring6[i] = pmem
        ring6[i].fitness = -i - 1
    for i in range(1000):
        deluge(ring6, ring6.occupied, inputs.normalise)
    assert ring6.num_filled == 1
    assert ring6[0].fitness == -1


def test_asteroid():
    inputs = Inputs()
    inputs.num_geoms = 0
    inputs.num_diheds = 3

    # smallest ring size is 5, so do testing on this thoroughly
    inputs.num_slots = 5
    ring3 = Ring(5, 0)
    for i in range(inputs.num_slots):
        pmem = Pmem(i, 0, 0, inputs.num_diheds)
        ring3[i] = pmem
    for i in range(1000):
        asteroid(ring3, ring3.occupied)
        assert ring3.num_filled >= 1
