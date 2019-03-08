"""Test the energy module of Kaplan."""

from kaplan.energy import run_energy_calc, MethodError, BasisSetError
from kaplan.inputs import Inputs
from kaplan.ring import Ring

def test_run_energy_calc():
    """Test the run_energy_calc function of the energy module."""
    # first make inputs
    inputs = Inputs()
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