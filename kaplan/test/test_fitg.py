"""Test the fitg module of Kaplan."""

from kaplan.fitg import run_energy_calc
from kaplan.inputs import Inputs
from kaplan.ring import Ring

def test_run_energy_calc():
    inputs = Inputs()
    inputs.update_inputs({
        "struct_type": "smiles",
        "struct_input": "CCC=CC=C",
        "num_mevs": 10,
        "num_geoms": 2,
        "num_slots": 10,
        "init_popsize": 5,
        "t_size": 2
    })
    ring = Ring(inputs.num_slots, inputs.init_popsize)
    print(ring[0])

test_run_energy_calc()