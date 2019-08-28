
from kaplan.output import run_output
from kaplan.pmem import Pmem
from kaplan.inputs import Inputs, read_input, InputError
from kaplan.ring import Ring
from kaplan.tournament import run_tournament
from kaplan.fitness import update_all_fitness

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
    print(inputs.obmol)
    ring = Ring(10, 5)
    print(ring)
    with assert_raises(InputError):
        # no fitness values have been set
        run_output(ring, 0, True)
    for i in range(1, inputs.num_mevs):
        update_all_fitness(ring)
        run_output(ring, i, True)
        run_tournament(ring, i)
    assert_raises(AssertionError, run_output, ring, "final", True)
    update_all_fitness(ring)
    run_output(ring, "last", True)
