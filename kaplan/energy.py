

# TODO: make these functions callable from a function
# that checks the program being used (i.e. psi4 vs horton
# vs gaussian)


import psi4
# link to relevant documentation
# http://www.psicode.org/psi4manual/1.2/psiapi.html

# make the output to the terminal be quiet
psi4.core.be_quiet()

# alternatively, can save output to a file
#psi4.core.set_output_file('random.txt', False)

# set RAM usage
# note: not sure this is called if done outside of
# the function call
RAM = "2 GB"
psi4.set_memory(RAM)


def run_energy_calc(geom, method="scf",basis="aug-cc-pVTZ",
                    restricted=False):
    """Run an energy calculation using psi4.

    Parameters
    ----------
    geom : str
        The string representing the molecular
        geometry for which to evaluate the energy.
        Note: this string should be generated using
        the prep_psi4_geom() function.

    method : str
        The quantum mechanical method to use.

    basis : str
        The basis set to use for the calculation.

    restricted : bool
        Default is False (runs an unrestricted
        calculation). If set to True, runs a
        restricted calculation.

    Raises
    ------
    AssertionError
        The method, basis, or geom variable
        is not a string.

    Returns
    -------
    A floating point number representing the energy
    of the given molecule, for the given basis set
    and method, in hartrees.

    Notes
    -----
    Restricted might not work for non-hf methods.

    """
    assert isinstance(method, str)
    assert isinstance(basis, str)
    assert isinstance(geom, str)
    molecule = psi4.geoemtry(geom)
    if restricted:
        psi4.set_options({"reference": "uhf"})
    energy = psi4.energy(method+'/'+basis)
    return energy

def prep_psi4_geom(coords, charge, multip):
    """Make a psi4 compliant geometry string.

    Parameters
    ----------
    coords : list(list)
        Atomic cartesian coordinates and atom types.
        Example for H_2:
        [['H', 0.0, 0.0, 0.0], ['H', 0.0, 0.0, 1.0]]
    charge : int
        The charge of the molecule.
    multip : int
        The multiplicity of the molecule.

    Returns
    -------
    A string as per psi4 input.

    """
    psi4_str = f"{charge} {multip}\n"
    for atom in coords:
        psi4_str += f"{atom[0]} {atom[1]} {atom[2]} {atom[3]}\n"
    return psi4_str
    
    


# jen messing around

#print('the thing')


#water = psi4.geometry("""
#O
#H 1 0.96
#H 1 0.96 2 104.5
#""")

#result = psi4.energy('scf/cc-pvdz')

#print(result)
#print(type(result))
#print(result.__dir__())


## another test

#psi4.set_memory("2 GB")

#methane = psi4.geometry("""
#C  0.0000    0.0000    0.0000
#H  0.5541    0.7996    0.4965
#H  0.6833   -0.8134   -0.2536
#H -0.7782   -0.3735    0.6692
#H -0.4593    0.3874   -0.9121
#""")

#result2 = psi4.energy('scf/cc-pvdz')
#print(result)
