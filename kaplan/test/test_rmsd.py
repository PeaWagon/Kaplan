
import os

from kaplan.rmsd import calc_rmsd
from vetee.xyz import Xyz

# directory for this test file
test_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'testfiles')

def test_calc_rmsd():
    # test the A1 and A3 hydrogen molecules
    mol1 = Xyz(os.path.join(test_dir, "H2-1A.xyz"))
    mol3 = Xyz(os.path.join(test_dir, "H2-3A.xyz"))
    assert calc_rmsd(mol1.coords, mol3.coords) == 1.0
    # test translated/rotated hydrogen
    mol1tr = Xyz(os.path.join(test_dir, "H2-1A-transrot.xyz"))
    assert calc_rmsd(mol1.coords, mol1tr.coords) == 0.0
    # test same molecule twice
    assert calc_rmsd(mol1.coords, mol1.coords) == 0.0

if __name__ == "__main__":
    test_calc_rmsd()
