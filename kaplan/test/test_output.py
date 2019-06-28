
from kaplan.output import energy_barplot, energy_rmsd_scatter, dihedrals_heatmap
from kaplan.pmem import Pmem
from kaplan.inputs import Inputs, read_input
from kaplan.ring import Ring

def test_energy_rmsd_scatter():
    inputs = Inputs()
    inputs.update_inputs({
        "struct_input": "propane"
    })
    print(inputs.num_dihed)
    p1 = Pmem(0, 0, 5, 3)
    p1.set_fitness()
    p2 = Pmem(0, 0, 5, 3)
    p3 = Pmem(0, 0, 5, 3)
    result = energy_barplot(p1)
    print(result)

def test_energy_barplot():
    pass


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
    print(inputs.min_diheds)

test_dihedrals_heatmap()
#test_read_output()

#test_energy_barplot()
#test_energy_rmsd_scatter()