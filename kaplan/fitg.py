

from kaplan.energy import 
from kaplan.rmsd import 

def get_fitness(pmem, fit_form, coef_energy, coef_rmsd):
    pass

def calc_energy(parser):
    # TODO: add parser attribute prog
    # so we can do this:
    # if parser.prog != 'psi4': blah blah
    if parser.prog != 'psi4':
        raise NotImplementedError("Only psi4 is supported at this time.")
    input_geom = prep_psi4_geom(parser)
    run_energy_calc(input_geom, parser.method, parser.basis)



def sum_energies(



