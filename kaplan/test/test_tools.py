"""Test the tools module of kaplan."""

import os

from numpy.testing import assert_raises

from kaplan.inputs import Inputs, InputError
from kaplan.pmem import Pmem
from kaplan.ring import Ring
from kaplan.tools import TEST_DIR, profile_function, analyse_profile,\
    plot_2d, make_2d, energy_barplot, energy_rmsd_scatter, dihedrals_heatmap
from kaplan.fitness import set_absolute_fitness


# see what amino acids look like since I will be using
# these structures frequently
names = [
    "asparagine", "glutamine", "aspartate", "glycine",
    "tryptophan", "cysteine", "threonine", "alanine",
    "isoleucine", "leucine", "tyrosine", "glutamate",
    "proline", "histidine", "lysine", "serine",
    "arginine", "valine", "methionine", "phenylalanine",
]


def test_profile_function():
    """Test the profile_function function from kaplan.tools.

    Notes
    -----
    This test shows the implementation for profiling a
    function. Just have to give:
    (1) Name of function to profile.
    (2) Name of dump file.
    (3) arguments followed by keyword arguments
    Use analyse_profile to convert the pstats object
    to a plain text file.

    """
    def my_adder(a, b):
        return a + b

    def my_list_adder(lista, listb):
        total = 0
        for i, j in zip(lista, listb):
            total += my_adder(i, j)
        return total

    profile_function(my_list_adder, "test_pf.dmp", range(10_000), range(10_000))
    analyse_profile("test_pf.dmp", "test_pf.log")
    # clean the profile test
    os.remove("test_pf.dmp")
    os.remove("test_pf.log")


def test_plot2d():
    inputs = Inputs()
    inputfile = os.path.join(TEST_DIR, "caffeine.xyz")
    inputs.update_inputs({
        "struct_input": inputfile,
        "struct_type": "xyz",
        "charge": 0,
        "multip": 1,
    })
    plot_2d(inputfile)

    inputfile = os.path.join(TEST_DIR, "1,3-butadiene.xyz")
    inputs.update_inputs({
        "struct_input": inputfile,
        "struct_type": "xyz",
        "charge": 0,
        "multip": 1,
    })
    plot_2d(inputfile)

    for name in names:
        inputs.update_inputs({
            "struct_input": name,
        })
        plot_2d()


def test_make2d():
    inputs = Inputs()
    for name in names:
        inputs.update_inputs({
            "struct_input": name,
        })
        make_2d()


def test_energy_barplot():
    # mostly have to check the test output manually here
    # test barplot works for regular pmem
    inputs = Inputs()
    inputs.update_inputs({
        "struct_input": "propane",
        "num_geoms": 10,
        "normalise": False,
    })
    pmem = Pmem(0, 0, 10, inputs.num_diheds)
    pmem.setup(major=False)
    set_absolute_fitness(pmem, inputs.fit_form, inputs.coef_energy, inputs.coef_rmsd)
    energy_barplot(pmem)

    # test error is raised when pmem has no energies
    pmem2 = Pmem(4, 0, 10, inputs.num_diheds)
    assert_raises(InputError, energy_barplot, pmem2)
    pmem2.setup(major=False)
    set_absolute_fitness(pmem2, inputs.fit_form, inputs.coef_energy, inputs.coef_rmsd)
    pmem2.energies[3] = None
    pmem2.energies[4] = None
    pmem2.energies[0] = None
    # test barplot works when some energies are None and that the conf
    # numbers are correctly assigned
    # also make sure two separate plots are generated and that the names
    # correspond to the pmem locations
    energy_barplot(pmem2)

    # test energy range works even if the range is:
    # (1) neg to pos
    # (2) neg to neg (tested above)
    # (3) pos to pos
    # also here test that image is correctly overwritten
    pmem3 = Pmem(4, 0, 10, inputs.num_diheds)
    for i in range(10):
        pmem3.energies[i] = -i + 5
    # also test here what happens when in and out units are the same
    energy_barplot(pmem3, outunits="Ha")
    pmem4 = Pmem(6, 0, 10, inputs.num_diheds)
    for i in range(10):
        pmem4.energies[i] = i
    energy_barplot(pmem4)


def test_energy_rmsd_scatter():
    # tests must be checked manually here
    inputs = Inputs()
    inputs.update_inputs({
        "struct_input": "propane",
        "num_geoms": 10,
        "normalise": False,
        "prog": "openbabel",
    })
    # test sizes 1-15 to make sure scaling factor allows all values
    # to be shown; 1 should not work since you need 2 values to show delta/rmsd
    assert_raises(InputError, energy_rmsd_scatter, Pmem(0, 0, 1, inputs.num_diheds))
    for i in range(2, 16):
        pmem = Pmem(i, 0, i, inputs.num_diheds)
        pmem.setup(major=False)
        energy_rmsd_scatter(pmem)


def test_dihedrals_heatmap():
    inputs = Inputs()
    inputs.update_inputs({
        "struct_input": "cyclohexane",
        "num_slots": 50,
        "init_popsize": 50,
        "no_ring_dihed": False,
    })
    r = Ring(50, 50)
    print(r)
    heatmap = dihedrals_heatmap(r)
    print(heatmap)
    print(inputs.diheds)
