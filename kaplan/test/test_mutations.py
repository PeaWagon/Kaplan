"""Test the mutations module from Kaplan."""

from numpy.testing import assert_raises

from kaplan.mutations import generate_children

# num muts num swaps


def test_generate_children():
    """Test the generate_children function from the mutations module."""
    parent1 = [[-1, -2, -3, -4, -5], [-1, -2, -3, -4, -5], [-1, -2, -3, -4, -5]]
    parent2 = [[-6, -7, -8, -9, -10], [-6, -7, -8, -9, -10], [-6, -7, -8, -9, -10]]
    # no changes are applied
    child1, child2 = generate_children(parent1, parent2, 0, 0)
    assert child1 == parent1
    assert child2 == parent2
    # make maximum of one mutation (to each child, for each geom)
    child1, child2 = generate_children(parent1, parent2, 1, 0)
    # go through changes and assert maximum 6 changes were made
    num_changes = 0
    for i, geom in enumerate(child1):
        for j, dihedral in enumerate(geom):
            if dihedral != parent1[i][j]:
                assert 360 > dihedral >= 0
                num_changes += 1
            if child2[i][j] != parent2[i][j]:
                # make sure dihedral angle is valid
                assert 360 > child2[i][j] >= 0
                num_changes += 1
    assert num_changes <= 6
    # make maximum of one swap
    child1, child2 = generate_children(parent1, parent2, 0, 1)
    num_changes1 = 0
    for i, geom in enumerate(child1):
        if geom != parent1[i]:
            num_changes1 += 1
    num_changes2 = 0
    for i, geom in enumerate(child2):
        if geom != parent2[i]:
            num_changes2 += 1
    assert num_changes1 == num_changes2
    assert num_changes1 <= 1
    # AssertionError when num_muts > num_atoms
    # and when num_swaps > num_geoms
    # or if parents are not the same length
    with assert_raises(AssertionError):
        generate_children(parent1, [[1, 2, 3, 4, 5], [13, 4, 2, 2, 5]], 0, 0)
        generate_children(parent1, parent2, 10, 0)
        generate_children(parent1, parent2, 0, 100)
    # just run the function a few times
    generate_children(parent1, parent2, 5, 3)
    generate_children(parent1, parent2, 4, 3)
    generate_children(parent1, parent2, 5, 2)
    generate_children(parent1, parent2, 3, 2)
