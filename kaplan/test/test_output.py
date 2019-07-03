
from kaplan.output import run_output
from kaplan.pmem import Pmem
from kaplan.inputs import Inputs, read_input
from kaplan.ring import Ring
from kaplan.tournament import run_tournament

from numpy.testing import assert_raises


def test_run_output():
    inputs = Inputs()
    inputs.update_inputs({
        "struct_input": "cysteine",
        "num_mevs": 5,
        "prog": "openbabel",
        "num_geoms": 1,
        "num_slots": 10,
        "init_popsize": 5,
    })
    ring = Ring(10, 5)
    print(ring)
    for i in range(inputs.num_mevs):
        run_output(ring, i, True)
        run_tournament(ring, i)
    assert_raises(AssertionError, run_output, ring, "final", True)
    run_output(ring, "last", True)
