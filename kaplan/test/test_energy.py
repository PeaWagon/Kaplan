"""Test the energy module of Kaplan."""

from kaplan.energy import run_energy_calc, MethodError, BasisSetError
from kaplan.inputs import Inputs

def test_run_energy_calc():
    """Test the run_energy_calc function of the energy module."""
    # first make inputs
    inputs = Inputs()

