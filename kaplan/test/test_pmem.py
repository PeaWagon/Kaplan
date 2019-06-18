from kaplan.pmem import Pmem
from kaplan.inputs import Inputs

# profiling
import cProfile, pstats, io

def test_pmem():
    inputs = Inputs()
    inputs.update_inputs({"struct_input": "propane"})
    p = Pmem(0, 0, 5, 2)
    #print(p)
    #for conf in p:
    #    print(conf)
    for i, _ in enumerate(p):
        p.set_energy_get_coords(i)
    #p.set_energy_get_coords(0)
    assert all(p.energies[i] is not None for i in range(len(p.energies)))
    #print(p.energies)
    #print(p)
    #p.set_fitness()

def profile_pmem():
    pr = cProfile.Profile()
    pr.enable()
    test_pmem()
    pr.disable()
    s = io.StringIO()
    sortby = pstats.SortKey.CUMULATIVE
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats()
    print(s.getvalue())



profile_pmem()