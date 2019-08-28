from kaplan.pmem import Pmem
from kaplan.inputs import Inputs
from kaplan.fitness import set_absolute_fitness
from copy import deepcopy

# profiling
import cProfile
import pstats
import io


def test_pmem():
    inputs = Inputs()
    inputs.update_inputs({"struct_input": "propane"})
    p = Pmem(0, 0, 5, 2)
    print(p)
    for conf in p:
        print(conf)
    p.setup(major=False)
    assert all(p.energies[i] is not None for i in range(len(p.energies)))
    assert all(p.rmsds[i][2] is not None for i in range(p.num_pairs))
    set_absolute_fitness(p, inputs.fit_form, inputs.coef_energy, inputs.coef_rmsd)
    assert p.fitness is not None
    print(p)


def profile_pmem():
    pr = cProfile.Profile()
    pr.enable()
    test_pmem()
    pr.disable()
    s = io.StringIO()
    sortby = pstats.SortKey.CUMULATIVE
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats()
    print(s.getvalue()[:int(0.05 * len(s.getvalue()))])


def test_pmem_get_geometry():
    inputs = Inputs()
    inputs.update_inputs({
        "struct_input": "asparagine",
        "num_geoms": 5,
    })
    p1 = Pmem(0, 0, inputs.num_geoms, inputs.num_diheds)
    p2 = deepcopy(p1)
    p1.setup(major=False)
    p2.setup(major=False)
    for energy1, energy2 in zip(p1.energies, p2.energies):
        assert energy1 == energy2
