"""Test the energy module of Kaplan."""

import os
import pybel
import openbabel

import numpy as np

from kaplan.energy import run_energy_calc, MethodError,\
    BasisSetError, obabel_energy, prep_psi4_geom,\
    psi4_energy_calc
from kaplan.inputs import Inputs
from kaplan.ring import Ring
from kaplan.tools import TEST_DIR, amino_acids
from kaplan.geometry import get_coords


def test_run_energy_calc():
    """Test the run_energy_calc function of the energy module."""
    # inputs for E-Z isomers of 1-phenyl-1,3-butadiene
    inputs = Inputs()
    inputs.update_inputs({
        "struct_input": os.path.join(TEST_DIR, "trans-1-phenyl-1,3-butadiene.xyz"),
        "struct_type": "xyz",
        "charge": 0,
        "multip": 1,
        "init_popsize": 5,
        "num_slots": 20,
        "mating_rad": 3
    })
    Ring(20, 5)

    # first make inputs
    inputs.update_inputs({
        "struct_input": "CC=CCCl",
        "struct_type": "smiles",
        "num_geoms": 3,
        "num_mevs": 10,
        "num_slots": 20,
        "init_popsize": 2,
        "mating_rad": 2
    })
    run_energy_calc(inputs.coords)

    # try with openbabel as prog
    inputs.update_inputs({
        "struct_input": "butane",
        "prog": "openbabel",
    })
    run_energy_calc(inputs.coords)


def test_prep_psi4_geom():
    """Test the prep_psi4_geom function from the energy module."""
    inputs = Inputs()
    inputs.update_inputs({"struct_input": "propane"})
    geom = prep_psi4_geom(inputs.coords, inputs.atomic_nums, inputs.charge, inputs.multip)
    assert isinstance(geom, str)
    geom_as_list = geom.strip().split("\n")
    assert geom_as_list[0] == f"{inputs.charge} {inputs.multip}"
    assert geom_as_list[-1] == "units angstrom"
    for atom in range(1, 12):
        atom_as_list = geom_as_list[atom].split(" ")
        assert atom_as_list[0].isalpha()
        float(atom_as_list[1])
        float(atom_as_list[2])
        float(atom_as_list[3])


def test_psi4_energy_calc():
    """Test the psi4_energy_calc function from the energy module."""
    test_mols = ["propane", "butane", "isobutane", "glycine", "benzene"]
    inputs = Inputs()
    for mol in test_mols:
        inputs.update_inputs({
            "struct_input": mol, "prog": "psi4",
            "method": "hf", "basis": "STO-3G",
            "no_ring_dihed": False,
        })
        geom = prep_psi4_geom(inputs.coords, inputs.atomic_nums, inputs.charge, inputs.multip)
        energy = psi4_energy_calc(geom, inputs.method, inputs.basis, "4 GB")
        assert energy < 0
        print(f"{mol} = {energy} hartrees")


def test_obabel_energy():
    # test to make sure the energy calculation with
    # MMFF94 is the same as the terminal command
    # results via obenergy filename.sdf

    # print this and run in terminal in testfiles directory
    # from kaplan.tools import amino_acids
    # for aa in amino_acids:
    #     print(f"obenergy {aa}.sdf | grep 'TOTAL ENERGY'&& ", end=" ")
    terminal_energies = {
        "asparagine": -9.2819,
        "glutamine": -7.44254,
        "aspartate": -6.60361,
        "glycine": 23.01667,
        "tryptophan": 53.58786,
        "cysteine": 42.61119,
        "threonine": 54.84375,
        "alanine": 26.06462,
        "isoleucine": 29.37363,
        "leucine": 36.21058,
        "tyrosine": 46.26394,
        "glutamate": -2.60574,
        "proline": 29.27256,
        "histidine": 41.68975,
        "lysine": 21.18729,
        "serine": 46.82625,
        "arginine": -101.21649,
        "valine": 32.33713,
        "methionine": 29.84991,
        "phenylalanine": 56.51525,
    }
    for aa in amino_acids:
        sdf_file = os.path.join(TEST_DIR, f"{aa}.sdf")
        assert os.path.isfile(sdf_file)
        mol = pybel.readfile("sdf", sdf_file)
        mol = mol.__next__()
        obmol = mol.OBMol
        energy_from_func = obabel_energy("mmff94", obmol)
        assert np.allclose(energy_from_func, terminal_energies[aa])
