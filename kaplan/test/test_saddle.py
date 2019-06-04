
import os

from saddle import internal
from saddle.coordinate_types import BondLength

ANG_TO_AU = lambda x : x/0.52917721067
AU_TO_ANG = lambda x : x*0.52917721067

TEST_DIR = "testfiles/"

test1 = os.path.join(TEST_DIR, "3-methylnona-3,4-diene.xyz")
test2 = os.path.join(TEST_DIR, "2-pentanol.xyz")
test3 = os.path.join(TEST_DIR, "caffeine.xyz")

#mol = internal.Internal.from_file(test2)
#mol.auto_select_ic(minimum=True)
#for ic in mol.ic:
#    print(ic)

mol = internal.Internal.from_file(test3)
mol.auto_select_ic(minimum=True)
for i, coord in enumerate(mol.coordinates):
    print([str(AU_TO_ANG(x))[:6] for x in coord])
for ic in mol.ic:
    if isinstance(ic, BondLength):
        print(ic)

