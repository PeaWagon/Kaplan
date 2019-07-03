
import os

from saddle import internal
from saddle.coordinate_types import BondLength

from kaplan.tools import TEST_DIR
from kaplan.geometry import geometry_units


def test_internal():
    """Test the internal module from saddle/GOpt."""
    files = [
        "3-methylnona-3,4-diene.xyz",
        "2-pentanol.xyz",
        "caffeine.xyz",
    ]
    for f in files:
        test = os.path.join(TEST_DIR, f)
        mol = internal.Internal.from_file(test)
        mol.auto_select_ic(minimum=True)
        for ic in mol.ic:
            print(ic)
            if isinstance(ic, BondLength):
                print(ic)
        for coord in mol.coordinates:
            print([str(geometry_units["atomic units"]["angstroms"](x))[:6] for x in coord])
