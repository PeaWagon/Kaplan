"""
Here is a list of all of the functions and objects in kaplan
This list is imported when the user writes "from kaplan import *"
"""

import pkg_resources

# enables use of kaplan.__version__
__version__ = pkg_resources.get_distribution("kaplan").version

from kaplan.energy import DEFAULT_INPUTS, BasisSetError, MethodError,\
                          run_energy_calc, prep_psi4_geom, psi4_energy_calc
from kaplan.control import run_kaplan
from kaplan.geometry import MAX_VALUE, MIN_VALUE, ANG_TO_AU, RAD_TO_DEGREES,\
                            get_coords
from kaplan.mutations import generate_children, mutate, swap, crossover, single_parent_mutation
from kaplan.output import run_output
from kaplan.pmem import Pmem
from kaplan.ring import RingEmptyError, RingOverflowError, Ring
from kaplan.rmsd import calc_rmsd
from kaplan.tournament import run_tournament, select_pmems, select_parents
from kaplan.inputs import Inputs, InputError, read_input
from kaplan.tools import TEST_DIR, generate_data, profile_function, analyse_profile, plot_2d, constants, units
from kaplan.extinction import apply_extinction, asteroid, agathic, deluge, plague