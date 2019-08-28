from numpy import allclose

from kaplan.fitness import list_replace, update_all_fitness
from kaplan.ring import Ring
from kaplan.inputs import Inputs
from kaplan.control import run_kaplan


def test_list_replace():
    """Test the list_replace function from the fitness module."""
    # basic implementation
    full = [1, 2, 3, 4, 5, 6]
    old = [1, 2, 3]
    new = [8, 9, 10]
    list_replace(full, old, new)
    assert full == [8, 9, 10, 4, 5, 6]

    # if duplicates exist
    full = [1, 2, 1, 2, 3, 5, 6]
    old = [1, 2, 3]
    new = [8, 9, 10]
    list_replace(full, old, new)
    assert full == [8, 9, 1, 2, 10, 5, 6]

    # if old is shorter than new
    full = [1, 2, 1, 2, 3, 5, 6]
    old = [1, 2]
    new = [8, 9, 10]
    list_replace(full, old, new)
    assert full == [8, 9, 1, 2, 3, 5, 6, 10]

    # if new is shorter than old
    full = [1, 2, 1, 2, 3, 5, 6]
    old = [1, 2, 3]
    new = [8, 9]
    list_replace(full, old, new)
    assert full == [8, 9, 1, 2, 5, 6]

    # with floating point values
    full = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    old = [1.0, 2.0, 3.0]
    new = [8.0, 9.0, 10.0]
    list_replace(full, old, new)
    assert allclose(full, [8.0, 9.0, 10.0, 4.0, 5.0, 6.0])


def test_update_all_fitness():
    # want to make sure fitness calculation
    # works with smallest possible size
    inputs = Inputs()
    inputs.update_inputs({
        "struct_input": "ethane",
        "init_popsize": 1,
        "num_slots": 10,
        "num_geoms": 1,
        "normalise": True,
    })
    run_kaplan(inputs)
    inputs.update_inputs({
        "struct_input": "ethane",
        "init_popsize": 1,
        "num_slots": 10,
        "num_geoms": 1,
        "normalise": False,
    })
    run_kaplan(inputs)
    inputs.update_inputs({
        "struct_input": "ethane",
        "init_popsize": 1,
        "num_slots": 10,
        "num_geoms": 1,
        "normalise": True,
        "stop_at_conv": True,
    })
    run_kaplan(inputs)
    inputs.update_inputs({
        "struct_input": "ethane",
        "init_popsize": 1,
        "num_slots": 10,
        "num_geoms": 1,
        "normalise": False,
        "stop_at_conv": True,
    })
    run_kaplan(inputs)
    inputs.update_inputs({
        "struct_input": "butane",
        "init_popsize": 1,
        "num_slots": 10,
        "num_geoms": 1,
        "normalise": False,
        "stop_at_conv": True,
    })
    run_kaplan(inputs)
