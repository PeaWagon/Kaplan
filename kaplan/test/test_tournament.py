
from kaplan.tournament import run_tournament, select_pmems, select_parents
from kaplan.ring import Ring
from vetee.parser import Parser

def test_run_tournament():
#self, num_geoms, num_atoms, num_slots,
#                 pmem_dist, fit_form, coef_energy, coef_rmsd,
#                 parser)
    #t_size, num_muts, num_swaps, ring,
#                   current_mev)

    parser = Parser()
    ring = Ring(3, 10, 50, 2, 0, 0.5, 0.5, parser)
    ring.fill(20, 0)
    for i in range(20):
        ring.pmems[i].fitness = i+20
    print('helo')
    print([ring.pmems[i].fitness for i in range(20)])
    run_tournament(4, 3, 1, ring, 0)

def test_select_pmems():
    pass

def test_select_parents():
    pass

if __name__ == "__main__":
    test_run_tournament()
    test_select_pmems()
    test_select_parents()
