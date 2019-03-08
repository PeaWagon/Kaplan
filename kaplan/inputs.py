"""Default values for inputs as "borg"
class structure.

If the _options dictionary is changed
(by instatiating another of) Inputs
or DefaultInputs, the _options is updated
for all users.

This class is not responsible for checking
initial convergence or inputs related to the
energy calculations (i.e. it does not check
if the basis set and/or method are available
for the given program). These responsibilities
will be shifted to the energy module.

"""


import numpy as np
import vetee

from saddle.internal import Internal


# 1 "angstrom" = 1.8897261339213 "atomic unit (length)"
#AU = 1.8897261339213
# https://github.com/psi4/psi4/blob/fbb2ff444490bf6b43cb6e027637d8fd857adcee/psi4/include/psi4/physconst.h
# inverse of bohr to angstrom from psi4
AU = 1/0.52917721067


class InputError(Exception):
    """Raised when invalid input is given."""


class DefaultInputs:


    # programs available in kaplan
    _avail_progs = {"psi4"}

    # available structure inputs for vetee to use to generate coordinates
    _avail_structs = {"smiles", "com", "xyz", "glog", "name", "cid"}

    # available fitness function formats
    _avail_fit_form = [0]

    _options = {
        # mol inputs
        "prog": "psi4",
        "basis": "sto-3g",
        "method": "hf",
        "struct_input": None,
        "struct_type": None,
        "num_atoms": None,
        "charge": None,
        "multip": None,
        # ga inputs
        "num_mevs": 10_000,
        "num_slots": 100,
        "init_popsize": 10,
        "pmem_dist": 3,
        "t_size": 7,
        # TODO: num_geoms should be defaulted based on a function
        # i.e. the num_geoms is guessed by the program
        "num_geoms": 5,
        # num_swaps and num_muts can have defaults based
        # on the num_geoms selection/default and the num_atoms
        # default for num_swaps is equal to num_geoms when not set
        "num_swaps": None,
        # default for num_muts is equal to num_atoms/3 when not set
        "num_muts": None,
        "fit_form": 0,
        "coef_energy": 0.5,
        "coef_rmsd": 0.5,
        # geometry specification for GOpt
        # none of these should be set by the user
        "atomic_nums": None,
        "coords": None,
        # How many dihedrals can be used to represent
        # the molecule (assuming one dihedral per
        # rotatable bond).
        "num_dihed": None,
        # List of indices representing where the dihedrals
        # are in the list of internal coordinates for the
        # molecule.
        "dihed_indices": None,
        # parser object from vetee
        "parser": None,
        # extra is a dictionary that
        # contains keywords to use in the energy
        # calculations (such as convergence criteria)
        # these may be program specific
        "extra": {}
    }


    def __init__(self):
        self.__dict__ = self._options


class Inputs(DefaultInputs):


    def __init__(self):
        super().__init__()


    def _reset_to_defaults(self):
        """Resets the borg to its original state.
        
        Notes
        -----
        This function is meant for testing purposes
        only. It should not be called during program
        execution.
        
        """
        self.prog = "psi4"
        self.basis = "sto-3g"
        self.method = "hf"
        self.struct_input = None
        self.struct_type = None
        self.num_atoms = None
        self.charge = None
        self.multip = None
        self.num_mevs = 10_000
        self.num_slots = 100
        self.init_popsize = 10
        self.pmem_dist = 3
        self.t_size = 7
        self.num_geoms = 5
        self.num_swaps = None
        self.num_muts = None
        self.fit_form = 0
        self.coef_energy = 0.5
        self.coef_rmsd = 0.5
        self.parser = None
        self.atomic_nums = None
        self.coords = None
        self.num_dihed = None
        self.dihed_indices = None
        self.extra = {}

    def _check_input(self):
        """Check that the current inputs are valid.

        Raises
        ------
        InputError
            Can be raised due to improper type (string
            given when int was expected). Can be raised
            due to invalid option (prog does not exist,
            fit_form not available, etc.).

        Returns
        -------
        None

        """
        # these sets are key names that are expected
        # to be integers or floats respectively
        # keys should not appear in both sets
        expect_int = {
            "num_atoms",
            "charge",
            "multip",
            "num_mevs",
            "num_slots",
            "init_popsize",
            "pmem_dist",
            "t_size",
            "num_geoms",
            "num_swaps",
            "num_muts",
            "num_dihed"
        }
        expect_float = {"coef_energy", "coef_rmsd"}
        # make sure all values have been set
        # and are set to the correct type
        for arg, val in self._options.items():
            if val is None:
                # defaults that are dependent on molecule size
                # or other inputs
                # val has to be reassigned otherwise
                # a ValueError occurs when int(None) is performed
                if arg == "num_swaps":
                    self.num_swaps = self.num_geoms
                    val = self.num_swaps
                elif arg == "num_muts":
                    self.num_muts = int(self.num_dihed/2)
                    val = self.num_muts
                else:
                    raise InputError(f"Missing required input argument: {arg}")
            # int(float) doesn't throw an error, but int("float") does
            # want to only accept integers and whole-number floats
            if arg in expect_int:
                if isinstance(val, float) and val % 1 != 0.0:
                    raise InputError(f"Expected integer argument: {arg} = {val}")
                try:
                    setattr(self, arg, int(val))
                except ValueError:
                    raise InputError(f"Could not convert to int: {arg} = {val}")
            elif arg in expect_float:
                try:
                    setattr(self, arg, float(val))
                except ValueError:
                    raise InputError(f"Could not convert to float: {arg} = {val}")

        # update parser charge and multip if not already set
        if self.parser.charge is None:
            self.parser.charge = self.charge
        if self.parser.multip is None:
            self.parser.multip = self.multip

        # only program currently available is psi4
        # if a program is added, add it to _avail_progs list in
        # DefaultInputs class
        assert self.prog in self._avail_progs
        # fit form can only be 0
        # if new fit forms are added, add them to _avail_fit_form
        # list in DefaultInputs class
        assert self.fit_form in self._avail_fit_form
        assert self.num_atoms > 3
        # multiplicity cannot be negative 2S+1 (where S is spin)
        assert self.multip > 0
        assert self.init_popsize > 0
        assert 0 < self.num_slots >= self.init_popsize
        assert self.num_mevs > 0
        assert self.num_geoms > 0
        assert 0 <= self.num_swaps <= self.num_geoms
        assert 0 <= self.num_muts <= self.num_dihed
        assert 0 <= self.pmem_dist < self.num_slots/2
        assert self.coef_energy >= 0
        assert self.coef_rmsd >= 0
        # t_size must be at least 2 (for 2 parents)
        assert 2 <= self.t_size <= self.init_popsize


    def _update_parser(self):
        """Updates the parser object (from vetee).

        Notes
        -----
        The following are required to be set:
            1. struct_type
            2. struct_input
            3. method
            4. basis
        The following are set to vetee defaults
        if left as None:
            1. charge
            2. multip
            3. num_atoms

        """
        try:
            if self.struct_type is None:
                raise InputError()
            assert self.struct_type in self._avail_structs
            # check the structure file exists (if applicable)
            if self.struct_type in ("xyz", "glog", "com"):
                with open(self.struct_input, "r"):
                    pass
            # make parser object using vetee
            if self.struct_type == "xyz":
                self.parser = vetee.xyz.Xyz(self.struct_input)
            elif self.struct_type == "com":
                self.parser = vetee.com.Com(self.struct_input)
            elif self.struct_type == "glog":
                self.parser = vetee.glog.Glog(self.struct_input)
            elif self.struct_type in ("smiles", "cid", "name"):
                self.parser = vetee.structure.Structure(self.struct_type,
                                                           self.struct_input)
        except InputError:
            raise InputError(f"Missing required input argument: struct_type")
        except AssertionError:
            raise InputError(f"Invalid structure type: {self.struct_type}")
        except FileNotFoundError:
            raise FileNotFoundError(f"No such struct_input file: {self.struct_input}")
        except ValueError:
            raise InputError(f"\n---Error generating parser object---\
                               \nStructure type: {self.struct_type}\
                               \nStructure input: {self.struct_input}")

        # update basis set and method
        self.parser.parse_gkeywords(f"#{self.method} {self.basis}")
        # if the user did not input charge/multip/num_atoms,
        # then try to use the defaults from vetee
        # if any of these values is set to None,
        # an input error will be raised by the update_inputs function
        if self.charge is None:
            self.charge = self.parser.charge
        if self.multip is None:
            self.multip = self.parser.multip
        if self.num_atoms is None:
            self.num_atoms = len(self.parser.coords)
        # make sure they agree for number of atoms
        else:
            assert self.num_atoms == len(self.parser.coords)


    def _update_coords(self):
        """Update the coords and atomic_nums attributes.
        
        Notes
        -----
        * coords must be in atomic units. It is an (n,3)
        numpy array, where n is the number of atoms.
        * atomic_nums is a list representing
        the atomic numbers of the atoms in the same
        order that they appear in the coords array.
        
        """
        self.coords = np.zeros((self.num_atoms, 3), float)
        for i, atom in enumerate(self.parser.coords):
            self.coords[i][0] = AU * atom[1]
            self.coords[i][1] = AU * atom[2]
            self.coords[i][2] = AU * atom[3]
        self.atomic_nums = np.array(vetee.periodic_table(\
            [self.parser.coords[i][0] for i in range(self.num_atoms)]))
        
        # now determine number of dihedral angles and where they
        # are in the internal coordinates array
        self.num_dihed = 0
        self.dihed_indices = []
        mol = Internal(self.coords, self.atomic_nums,
                       self.charge, self.multip)
        # generate internal coordinates, with one
        # dihedral per rotatable bond
        mol.auto_select_ic(minimum=True)
        for i, intern_coord in enumerate(mol.ic):
            # if the internal coordinate involves 4 atoms,
            # it is a dihedral angle
            if intern_coord._coordinates.shape == (4,3):
                self.num_dihed += 1
                self.dihed_indices.append(i)
        
        # check that number of dihedral angles is not 0
        if self.num_dihed == 0:
            raise InputError("No dihedral angles could be generated for the input molecule.")


    def update_inputs(self, input_dict):
        """Update the inputs from their default values.

        Parameters
        ----------
        input_dict : dict
            Key-value pairs that should be changed
            in the self._options dictionary.

        Notes
        -----
        This function should only be called once during
        execution of Kaplan.

        Returns
        -------
        None

        """
        # for safety, set the Inputs options
        # to their defaults every time this is called
        self._reset_to_defaults()
        # make sure it is actually a dictionary
        # make sure dict is non-empty
        assert isinstance(input_dict, dict) and len(input_dict) > 0
        # iterate over and check inputs
        for arg, val in input_dict.items():
            # arg can be case insensitive
            arg = arg.lower()
            # make sure no new keys are added
            if arg not in self._options:
                raise InputError(f"Invalid input: {arg} = {val}")
            # make sure lowercase is used
            # unless it's for structure (i.e.
            # SMILES strings are case sensitive)
            if arg != "struct_input":
                try:
                    val = val.lower()
                # value was integer or float
                except AttributeError:
                    pass
            # overwrite default value
            setattr(self, arg, val)

        # now update parser object
        self._update_parser()

        # update coordinates using parser object
        self._update_coords()

        # check that the inputs are valid
        self._check_input()
