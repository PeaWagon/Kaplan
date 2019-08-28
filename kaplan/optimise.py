"""

The purpose of this module is to optimise coordinates
locally using numerical methods. Currently, Openbabel
is used to apply conjugate gradients to an Openbabel
molecule object.

In future versions, quantum chemistry methods such
as self-consistent field can be implemented to optimise
local coordinates. Possible implementation is available
here:
http://www.psicode.org/psi4manual/master/optking.html

"""


import pybel
import psi4

from psi4.driver import constants
from openbabel import OBFF_LOGLVL_LOW

from numpy import allclose

from kaplan.inputs import Inputs, hardware_inputs
from kaplan.energy import MethodError, BasisSetError,\
    prep_psi4_geom, run_energy_calc
from kaplan.geometry import get_coords, set_coords


def optimise_coords(coords, major):
    """Locally optimise a set of coordinates.

    Parameters
    ----------
    coords : np.array((num_atoms, 3), float)
        The coordinates to optimise. The atomic
        numbers and other items are taken from
        the current inputs.
    major : bool
        False means to do a minor
        optimisation (called during evolution). If
        True, does a major optimisation (designed for
        initial coordinates and final conformer
        geometries). Note that minor geometry
        optimisation is done with Openbabel during
        psi4 jobs - only the energy is evalutated
        using pis4.

    Returns
    -------
    tuple(float, np.array((num_atoms, 3), float)
        The energy and coordinates for the optimised
        structure.

    """
    assert isinstance(major, bool)
    inputs = Inputs()
    if not major and inputs.prog == "openbabel":
        # give obmol the geometry to be optimised
        set_coords(inputs.obmol, coords)
        energy, opt_coords = obabel_geometry_opt(
            inputs.method, inputs.obmol, inputs.minor_tolerance,
            inputs.minor_maxsteps, inputs.minor_sampling
        )
    elif not major and inputs.prog == "psi4":
        set_coords(inputs.obmol, coords)
        _, opt_coords = obabel_geometry_opt(
            "mmff94", inputs.obmol, inputs.minor_tolerance,
            inputs.minor_maxsteps, inputs.minor_sampling
        )
        energy = run_energy_calc(opt_coords)
    elif major and inputs.prog == "openbabel":
        set_coords(inputs.obmol, coords)
        energy, opt_coords = obabel_geometry_opt(
            inputs.method, inputs.obmol, inputs.major_tolerance,
            inputs.major_maxsteps, inputs.major_sampling
        )
    elif major and inputs.prog == "psi4":
        geom_str = prep_psi4_geom(
            coords, inputs.atomic_nums, inputs.charge, inputs.multip
        )
        energy, opt_coords = psi4_geometry_opt(
            geom_str, inputs.method, inputs.basis, hardware_inputs["psi4"]["RAM"]
        )
    else:
        raise NotImplementedError(f"Program {inputs.prog} not found.")
    return energy, opt_coords


def psi4_geometry_opt(geom_str, method, basis, ram):
    """Optimise a geometry using psi4.

    Parameters
    ----------
    geom_str : str
        The geometry string with which to construct
        the molecule being optimised. Should be made
        from prep_psi4_geom function from the kaplan
        energy module.
    method : str
        The method to use for optimisation.
    basis : str
        The basis set to use for optimisation.
    ram : str
        The amount of random-access memory allowed
        for the computation. For example: "4 GB".

    Returns
    -------
    tuple(float, np.array((num_atoms, 3), float))
        Energy (in hartrees) and the optimised
        coordinates.

    """
    assert isinstance(method, str)
    assert isinstance(basis, str)
    psi4.core.be_quiet()
    psi4.set_memory(ram)
    # geometry inputs are angstroms
    psi4.set_options({"cfour_units": "angstrom"})
    mol = psi4.geometry(geom_str)
    try:
        energy = psi4.optimize(
            f"{method}/{basis}", return_wfn=False, molecule=mol
        )
        psi4_matrix = mol.geometry()  # units are bohr here
        psi4_matrix.scale(constants.get("bohr2angstroms"))
        coords = psi4_matrix.to_array()  # this is a numpy array
        psi4.core.clean()
        return energy, coords
    except psi4.driver.p4util.exceptions.ValidationError:
        raise MethodError(f"Invalid method for Psi4: {method}")
    except psi4.driver.qcdb.exceptions.BasisSetNotFound:
        raise BasisSetError(f"Invalid basis set for Psi4: {basis}")


def obabel_geometry_opt(ff, obmol, tolerance=1e-6, max_steps=2500,
                        sample_every=100, logging=False):
    """Optimise an openbabel molecule with a forcefield.

    Parameters
    ----------
    ff : str
        The forcefield to use. Should be in pybel.forcefields.
    obmol : openbabel OBMol object
        The openbabel molecule to be optimised.
    tolerance : float
        After each iteration of conjugate gradients, the
        energy is compared with the previous iteration. If
        this energy difference is less than the tolerance,
        then the energy and coordinates at that step are returned.
        If tolerance is negative, then all of the max_steps are run.
    max_steps : int
        How many iterations to perform if the tolerance
        condition is not satisfied.
    sample_every : int
        If tolerance is not negative, indicates how many
        iterations to perform before checking the energy
        for the tolerance.
    logging : bool
        Defaults to False, which means that no logging
        is written to the terminal. True means that the
        openbabel logging level LOW is switched on.

    Notes
    -----
    This function uses conjugate gradients, but
    steepest descent is also available from Openbabel
    (not used here). To switch, just change the code
    ConjugateGradients to SteepestDescent.

    Returns
    -------
    tuple(float, np.array((num_atoms, 3), float))
        Representing the final energy and the related
        coordinates (for the obmol).

    """
    # in case forcefield is uppercase
    ff = ff.lower()
    if ff not in pybel.forcefields:
        raise MethodError(f"Forcefield unavailable in Openbabel: {ff}")
    ff_instance = pybel.ob.OBForceField.FindForceField(ff)
    # this makes sure the pointer is valid
    assert ff_instance

    if logging:
        # can also change log level to HIGH for more output
        ff_instance.SetLogLevel(OBFF_LOGLVL_LOW)
        ff_instance.SetLogToStdErr()

    result = ff_instance.Setup(obmol)
    if not result:
        raise MethodError("Unable to setup forcefield for molecule")

    # if tolerance is negative, just run the optimiser with max steps
    if tolerance < 0:
        ff_instance.ConjugateGradients(max_steps)
        energy = ff_instance.Energy(False)
    # otherwise run, checking delta energies
    else:
        prev_energy = ff_instance.Energy(False)
        for _ in range(1, max_steps + 1, sample_every):
            ff_instance.ConjugateGradients(sample_every)
            energy = ff_instance.Energy(False)
            if allclose(energy, prev_energy, atol=tolerance, rtol=0):
                # convergence achieved within max_steps
                break
            else:
                prev_energy = energy

    # this actually updates the obmol coordinates with
    # those that have been optimised (otherwise no changes
    # are made to the object)
    ff_instance.GetCoordinates(obmol)
    opt_coords = get_coords(obmol)
    return energy, opt_coords
