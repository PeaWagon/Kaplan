"""Test the energy module of Kaplan."""

import os

from kaplan.energy import run_energy_calc, MethodError, BasisSetError
from kaplan.inputs import Inputs
from kaplan.ring import Ring
from kaplan.tools import TEST_DIR


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
