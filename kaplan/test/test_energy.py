"""Test the energy module of Kaplan."""

import os

from kaplan.energy import run_energy_calc, MethodError, BasisSetError
from kaplan.inputs import Inputs
from kaplan.ring import Ring


TEST_DIR = os.path.join(os.path.dirname(os.path.realpath("__file__")), 'testfiles')


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
        "t_size": 3
    })
    r = Ring(20, 5)
    
    # first make inputs
    inputs.update_inputs({
        "struct_input": "CC=CCCl",
        "struct_type": "smiles",
        "num_geoms": 3,
        "num_mevs": 10,
        "num_slots": 20,
        "init_popsize": 2,
        "t_size": 2
    })
    result = run_energy_calc(inputs.coords)
    print(result)

test_run_energy_calc()