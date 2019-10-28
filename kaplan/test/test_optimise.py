"""Test the optimise module of Kaplan."""

import os
import pybel
import numpy as np
from numpy.testing import assert_raises
from qcelemental.exceptions import ValidationError

from kaplan.tools import amino_acids, TEST_DIR
from kaplan.energy import obabel_energy, psi4_energy_calc,\
    prep_psi4_geom, run_energy_calc
from kaplan.optimise import optimise_coords, obabel_geometry_opt,\
    psi4_geometry_opt, opt_with_rdkit
from kaplan.geometry import get_coords
from kaplan.inputs import Inputs


def test_optimise_coords():
    """Test the optimise_coords function from the optimise module."""
    test_mols = ["aspartate", "isopropyl alcohol", "glycine"]
    inputs = Inputs()
    # do psi4 tests
    for mol in test_mols:
        inputs.update_inputs({
            "struct_input": mol, "prog": "psi4",
            "method": "hf", "basis": "STO-3G",
        })
        initial_energy = run_energy_calc(inputs.coords)
        energy, coords = optimise_coords(inputs.coords, False)
        assert energy < initial_energy
        assert not np.allclose(coords, inputs.coords)
        energy2, coords2 = optimise_coords(inputs.coords, True)
        assert energy2 < energy < initial_energy
        assert not np.allclose(coords, coords2)
        assert not np.allclose(coords2, inputs.coords)
    # do openbabel tests for major and minor
    for mol in test_mols:
        inputs.update_inputs({
            "struct_input": mol, "prog": "openbabel",
        })
        initial_energy = run_energy_calc(inputs.coords)
        energy_major, coords_major = optimise_coords(inputs.coords, True)
        energy_minor, coords_minor = optimise_coords(inputs.coords, False)
        assert energy_major < energy_minor < initial_energy
        assert not np.allclose(coords_major, coords_minor)
        assert not np.allclose(coords_major, inputs.coords)
        assert not np.allclose(coords_minor, inputs.coords)


def test_psi4_geometry_opt():
    """Test the psi4_geometry_opt function from the optimise module."""
    test_mols = ["propane", "isopropyl alcohol", "glycine"]
    inputs = Inputs()
    for mol in test_mols:
        inputs.update_inputs({
            "struct_input": mol, "prog": "psi4",
            "method": "hf", "basis": "STO-3G",
        })
        geom = prep_psi4_geom(inputs.coords, inputs.atomic_nums, inputs.charge, inputs.multip)
        energy = psi4_energy_calc(geom, inputs.method, inputs.basis, "4 GB")
        assert energy < 0
        print(f"{mol} = {energy} hartrees")
        energy2, coords = psi4_geometry_opt(geom, inputs.method, inputs.basis, "4 GB")
        print(f"{mol} = {energy2} hartrees")
        assert energy2 < energy
        # check that energy agrees when rerun with psi4_energy_calc
        geom2 = prep_psi4_geom(coords, inputs.atomic_nums, inputs.charge, inputs.multip)
        energy3 = psi4_energy_calc(geom2, inputs.method, inputs.basis, "4 GB")
        assert np.allclose(energy3, energy2)

    # note: you can use the molecule method save_string_xyz
    # get the coordinates of a psi4 molecule as a string
    # in xyz format
    # print(mol.save_string_xyz())
    # where mol is constructed via:
    # mol = psi4.geometry(geom_str)

    # see what happens with garbage coordinates
    bad_geom = """
    0 1
    C 0.0 0.0 0.0
    H 0.0 0.0 0.0
    Mg 0.0 0.0 0.0
    Au 0.1 0.3 0.4
    """
    # qcelemental.exceptions.ValidationError: Following atoms are too close:
    # [(1, 0, 0.0), (2, 0, 0.0), (2, 1, 0.0)]
    with assert_raises(ValidationError):
        psi4_geometry_opt(bad_geom, "hf", "sto-3g", "4 GB")


def test_obabel_geometry_opt():
    """Test the obabel_geometry_opt function from the energy module.

    Notes
    -----
    The obminimize tool does not stop
    iterations even when the energy has stabilised
    (within the specified crit value). This was
    reported already in Issue #1366 (github):
    https://github.com/openbabel/openbabel/issues/1366
    Furthermore, the energy of the Openbabel molecule
    does not decrease when read as a SMILES string
    and subsequently optimised. This behaviour is
    not reproducible when data is read from a file.
    This problem may have something to do with how
    many decimal places (from atomic coordinates) are
    considered when calculating energies versus
    generating geometries from SMILES.

    The results from running max steps (with tolerance
    set to None) are not equivalent to running with
    a numerical tolerance for the same molecule. If
    there are flat areas in the energy profile, then
    the optimiser may return without finding the
    global minimum.

    """
    # generated using following shell command:
    # obminimize -cg -n 2500 -ff MMFF94 -o sdf
    # $infile >> output.sdf 2>> $logfile
    # where $infile is the name of the amino
    # acid sdf file
    obminimise_energies = {
        "asparagine": -12.89555,
        "glutamine": -10.23275,
        "aspartate": -15.92273,
        "glycine": 20.9297,
        "tryptophan": 49.44199,
        "cysteine": 36.39594,
        "threonine": 50.72826,
        "alanine": 22.85702,
        "isoleucine": 26.42489,
        "leucine": 30.33185,
        "tyrosine": 41.85679,
        "glutamate": -7.60834,
        "proline": 24.84185,
        "histidine": 35.56765,
        "lysine": 18.75389,
        "serine": 43.93566,
        "arginine": -112.90437,
        "valine": 29.11269,
        "methionine": 23.3587,
        "phenylalanine": 52.35171,
    }
    for aa in amino_acids:
        sdf_file = os.path.join(TEST_DIR, f"{aa}.sdf")
        assert os.path.isfile(sdf_file)
        mol = pybel.readfile("sdf", sdf_file)
        mol = mol.__next__()
        obmol = mol.OBMol
        print()
        energy_from_func, new_coords = obabel_geometry_opt(
            "mmff94", obmol, tolerance=-1, max_steps=2500
        )
        print("from function:", energy_from_func)
        print("original terminal:", obminimise_energies[aa])
        assert np.allclose(new_coords, get_coords(obmol))
        assert np.allclose(obabel_energy("mmff94", obmol), energy_from_func)
        print(aa)
        print(f"energy diff: {energy_from_func - obminimise_energies[aa]}")
        print()
        assert np.allclose(energy_from_func, obminimise_energies[aa], atol=0.0001)


def test_opt_with_rdkit():
    """Test the opt_with_rdkit function from the optimise module."""
    inputs = Inputs()
    inputs.update_inputs({
        "struct_input": "butane",
    })
    # test that doing zero steps of optimisation returns the
    # same coordinates
    coords, result = opt_with_rdkit(inputs.obmol, 0)
    assert np.allclose(coords, inputs.coords)
    assert result == 1

    # test that optimisation converges for butane and changes coordinates
    opt_coords, opt_result = opt_with_rdkit(inputs.obmol, 500)
    assert not np.allclose(inputs.coords, opt_coords)
    assert opt_result == 0
