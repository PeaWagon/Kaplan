"""
This test file is for testing the ring iterator
and __str__ for pmem and ring objects.

"""

from kaplan.inputs import Inputs
from kaplan.ring import Ring
from kaplan.pmem import Pmem

def test_ring_iteration():
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
