"""This module uses psi4 to run energy calculations
for a given geometry."""

import os
import psi4

from kaplan.inputs import Inputs

# these values are used in energy calculations when no
# options are provided
# if a new program is added, add defaults for it
# how much RAM to use for psi4 calculations
# should be less than what your computer has available
DEFAULT_INPUTS = {"psi4": {"RAM": "4 GB"},
                  "other": {}
}


class BasisSetError(Exception):
    """Raised if basis set cannot be found for a program."""


class MethodError(Exception):
    """Raised if quantum chemical method is unavailable."""


def run_energy_calc(coords):
    """Setup and execute an energy calculation.

    Parameters
    ----------
    coords : list
        Has the format [[a1, x1, y1, z1],
        [a2, x2, y2, z2], ..., [an, xn, yn, zn]]
        Where a's are atom labels (i.e. H for
        hydrogen) and x, y, and z coordinates
        are given for n atoms. a's should be
        strings and xyz should be floats.

    Notes
    -----
    Other options can be added and incorporated.
        options : dict
        Requires the following:
        options["charge"] - molecule charge
        options["multip"] - molecule multiplicity
        These options can be set:
        options["prog"] - program (default is psi4)
        options["basis"] - basis set (default is sto-3g)
        options["method"] - quantum chem method (default is hf)
        options["RAM"] - memory to use (default is 4 GB)

    Returns
    -------
    Energy result as a floating point number.

    """
    inputs = Inputs()
    # if you added extra options, they are included here
    extras = {}
    for arg, val in inputs.extra.items():
        extras[arg] = val

    defaults = DEFAULT_INPUTS[inputs.prog]
    for key in defaults:
        if key not in extras:
            extras[key] = defaults[key]
    if inputs.prog == "psi4":
        geom_str = prep_psi4_geom(coords, inputs.charge, inputs.multip)
        energy = psi4_energy_calc(geom_str, inputs.method, inputs.basis,
                                  extras["RAM"])
        # energy is in hartrees here
        return energy
    raise NotImplementedError("Program not found.")


def psi4_energy_calc(geom, method, basis, ram, restricted=False):
    """Run an energy calculation using psi4.

    Parameters
    ----------
    geom : str
        The string representing the molecular
        geometry for which to evaluate the energy.
        Note: this string should be generated using
        the prep_psi4_geom() function.
    method : str
        The quantum mechanical method to use.
    basis : str
        The basis set to use for the calculation.
    ram : str
        Specifies how much memory to use in calculations.
        Should be something like "4 GB".
    restricted : bool
        Default is False (runs an unrestricted
        calculation). If set to True, runs a
        restricted calculation.

    Raises
    ------
    AssertionError
        The method or basis is not a string.
    BasisSetError
        The basis set did not work using psi4.
    MethodError
        The method did not work using psi4.

    Returns
    -------
    A floating point number representing the energy
    of the given molecule, for the given basis set
    and method, in hartrees.

    Notes
    -----
    Restricted might not work for non-hf methods.
    VERY untested.

    """
    psi4.core.be_quiet()
    psi4.set_memory(ram)
    assert isinstance(method, str)
    assert isinstance(basis, str)
    if restricted:
        psi4.set_options({"reference": "uhf"})
    try:
        return psi4.energy(method+'/'+basis, return_wfn=False)
    except psi4.driver.p4util.exceptions.ValidationError:
        raise MethodError(f"Invalid method: {method}")
    except psi4.driver.qcdb.exceptions.BasisSetNotFound:
        raise BasisSetError(f"Invalid basis set: {basis}")


def prep_psi4_geom(coords, charge, multip):
    """Make a psi4 compliant geometry string.

    Parameters
    ----------
    coords : list(list)
        Atomic cartesian coordinates and atom types.
        Example for H_2:
        [['H', 0.0, 0.0, 0.0], ['H', 0.0, 0.0, 1.0]]
    charge : int
        The charge of the molecule.
    multip : int
        The multiplicity of the molecule.

    Returns
    -------
    A string as per psi4 input.

    """
    psi4_str = f"\n{charge} {multip}\n"
    for atom in coords:
        psi4_str += f"{atom[0]} {atom[1]} {atom[2]} {atom[3]}\n"
    return psi4.geometry(psi4_str)
