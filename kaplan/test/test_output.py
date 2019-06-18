
from kaplan.output import energy_barplot, energy_rmsd_scatter
from kaplan.pmem import Pmem
from kaplan.inputs import Inputs, read_input

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

    

#test_read_output()

#test_energy_barplot()
#test_energy_rmsd_scatter()