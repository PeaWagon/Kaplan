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
import numpy as np

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
    # make a new array to hold the internal coordinates
    new_ic_values = mol.ic_values

    # keep track of the dihedral angle we are currently updating
    dihed_index = 0
    # update new dihedral angles
    for i, value in enumerate(new_ic_values):
        if i in inputs.dihed_indices:
            new_ic_values[i] = dihedrals[dihed_index]
            dihed_index += 1
    
    # debugging test check. can remove later
    assert dihed_index == len(inputs.dihed_indices)
    
    # set the target to these new values
    mol.set_target_ic(new_ic_values)
    mol.optimize_to_target_ic(method="BFGS", max_iter=500, hess=True)
    return mol.coordinates
