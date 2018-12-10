
import os

from vetee.xyz import Xyz
from numpy.testing import assert_raises

from kaplan.tournament import run_tournament, select_pmems, select_parents
from kaplan.ring import Ring, RingEmptyError

# directory for this test file
test_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'testfiles')

def test_run_tournament():
    # ring init:
    # num_geoms, num_atoms, num_slots, pmem_dist, fit_form,
    # coef_energy, coef_rmsd, parser
    mol = Xyz(os.path.join(test_dir, "1,3-butadiene.xyz"))
    mol.charge = 0
    mol.multip = 1
    ring = Ring(3, 10, 20, 2, 0, 0.5, 0.5, mol)
    ring.fill(3, 0)
    # tournament call:
    # t_size, num_muts, num_swaps, ring, current_mev
    # check that an empty ring error is raised when trying to
    # run a tournament size that is too big
    assert_raises(RingEmptyError, run_tournament, 7, 0, 0, ring, 0)
    run_tournament(2, 0, 0, ring, 0)

#    for i in range(20):
#        ring.pmems[i].fitness = i+20
#    print('helo')
#    print([ring.pmems[i].fitness for i in range(20)])
#    run_tournament(4, 3, 1, ring, 0)

def test_select_pmems():
    pass

def test_select_parents():
    pass

if __name__ == "__main__":
    test_run_tournament()
    test_select_pmems()
    test_select_parents()
