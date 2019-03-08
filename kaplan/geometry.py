"""This module is responsible for dealing with converting
geometries from z-matrix format to xyz format. It relies
on the external library, Vetee, to do most of the conversions,
along with the python wrapper for openbabel, pybel.

New changes to how dihedrals are used are being made.

Now this program will use saddle (GOpt) to generate
dihedral angles.

"""

import vetee
import openbabel
import pybel

from saddle.internal import Internal

from kaplan.inputs import Inputs


class GeometryError(Exception):
    """Error raises when one of the geometry functions fails."""


def get_geom_from_dihedrals(dihedrals):
    """Construct xyz coordinates using dihedral angles.

    Parameters
    ----------
    dihedrals : np.array(k)
        Where k is the number of dihedral angles
        in the molecule (as specified by GOpt).

    Returns
    -------
    coords : np.array(n, 3)
        xyz cartesian coordinates in atomic units.

    """
    inputs = Inputs()
    mol = Internal(inputs.coords, inputs.atomic_nums,
                   inputs.charge, inputs.multip)
    # generate internal coordinates, with one
    # dihedral per rotatable bond
    mol.auto_select_ic(minimum=True)
    # update new dihedral angles
    for i, index in enumerate(inputs.dihed_indices):
        mol.ic_values[index] = dihedrals[i]
    # set the target to these new values
    mol.set_target_ic(mol.ic_values)
    mol.optimize_to_target_ic()
    return mol.coordinates
