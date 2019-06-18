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


import os
import numpy as np

from vetee.coordinates import write_xyz

from kaplan.geometry import update_obmol, create_obmol,\
                            remove_ring_dihed, get_min_dihed, get_rings,\
                            get_struct_info, get_atomic_nums

from saddle.internal import Internal


class InputError(Exception):
    """Raised when invalid input is given."""


class DefaultInputs:


    # programs available in kaplan
    _avail_progs = {"psi4", "openbabel"}

    # available structure inputs for vetee to use to generate coordinates
    _avail_structs = {"smiles", "com", "xyz", "name", "cid", "inchi", "inchikey"}

    # available fitness function formats
    _avail_fit_form = [0]

    _options = {
        "output_dir": "pwd",    # where to store the output
        # mol inputs
        "prog": "psi4",         # program to use to run energy calculations
        "basis": "sto-3g",      # basis set to use in quantum calculations
        # if the prog is openbabel, use mmff94 as default forcefield
        # if the prog is psi4, use hf
        # if no method is specified, and com is input, then use method from com
        "method": None,         # method to run for energy evaluation
        "struct_input": None,   # value for input
        "struct_type": None,    # type of structural input (defaults to "name")
        "charge": None,         # total molecular charge
        "multip": None,         # molecule multiplicity
        "no_ring_dihed": True,  # remove ring dihedral angles
        # ga inputs
        "num_mevs": 100,
        "num_slots": 50,
        "init_popsize": 10,
        "mating_rad": 3,
        "num_geoms": 5,
        # default for num_swaps is equal to num_geoms//2 when not set
        "num_swaps": None,
        # default for num_muts is equal to (num_dihed*num_geoms)//5 when not set
        "num_muts": None,
        # default for num_cross is equal to num_geoms//2 when not set
        "num_cross": None,
        "fit_form": 0,
        "coef_energy": 0.5,
        "coef_rmsd": 0.5,
        # geometry specification for GOpt/openbabel
        # none of these should be set by the user
        "atomic_nums": None,    # atomic numbers by atom in the molecule
        # np.array((num_atoms,3), float) representing xyz coordinates
        # of the original input; preserved to ensure obmol retains
        # the same connections throughout optimisation
        "coords": None,         
        # How many dihedrals can be used to represent
        # the molecule (assuming one dihedral per
        # rotatable bond); dihedrals encased by rings
        # can be removed (set no_ring_dihed to True)
        "num_dihed": None,
        # list of tuples, where each tuple has 4 integers
        # representing atom indices that make up a minimum
        # dihedral angle
        "min_diheds": None,
        # extra is a dictionary that
        # contains keywords to use in the energy
        # calculations (such as convergence criteria)
        # these may be program specific
        "extra": {},
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
        self.output_dir = "pwd"
        self.prog = "psi4"
        self.basis = "sto-3g"
        self.method = None
        self.struct_input = None
        self.struct_type = None
        self.charge = None
        self.multip = None
        self.no_ring_dihed = True
        self.num_mevs = 100
        self.num_slots = 50
        self.init_popsize = 10
        self.mating_rad = 3
        self.num_geoms = 5
        self.num_swaps = None
        self.num_muts = None
        self.num_cross = None
        self.fit_form = 0
        self.coef_energy = 0.5
        self.coef_rmsd = 0.5
        self.atomic_nums = None
        self.coords = None
        self.num_dihed = None
        self.min_diheds = None
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
            "charge",
            "multip",
            "num_mevs",
            "num_slots",
            "init_popsize",
            "mating_rad",
            "num_geoms",
            "num_swaps",
            "num_muts",
            "num_cross",
            "num_dihed",
        }
        expect_float = {"coef_energy", "coef_rmsd"}
        expect_bool = {"no_ring_dihed"}
        # make sure all values have been set
        # and are set to the correct type
        for arg, val in self._options.items():
            if val is None:
                # defaults that are dependent on molecule size
                # or other inputs
                # val has to be reassigned otherwise
                # a ValueError occurs when int(None) is performed
                if arg == "num_swaps":
                    self.num_swaps = self.num_geoms//2
                    val = self.num_swaps
                elif arg == "num_cross":
                    self.num_cross = self.num_geoms//2
                    val = self.num_cross
                elif arg == "num_muts":
                    self.num_muts = (self.num_dihed*self.num_geoms)//5
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
            elif arg in expect_bool:
                try:
                    assert isinstance(val, bool)
                except ValueError:
                    raise InputError(f"Value should be True or False and not a string: {arg} = {val}")

        # only programs currently available are psi4 and openbabel
        # if a program is added, add it to _avail_progs list in
        # DefaultInputs class
        assert self.prog in self._avail_progs
        # fit form can only be 0
        # if new fit forms are added, add them to _avail_fit_form
        # list in DefaultInputs class
        assert self.fit_form in self._avail_fit_form
        # multiplicity cannot be negative 2S+1 (where S is spin)
        assert self.multip > 0
        assert self.init_popsize > 0
        assert 0 < self.num_slots >= self.init_popsize
        assert self.num_mevs > 0
        assert self.num_geoms > 0
        assert 0 <= self.num_swaps <= self.num_geoms
        assert 0 <= self.num_cross <= self.num_geoms
        assert 0 <= self.num_muts <= self.num_dihed*self.num_geoms
        assert 2 <= self.mating_rad <= self.num_slots//2
        assert self.coef_energy >= 0
        assert self.coef_rmsd >= 0


    def _update_geometry(self):
        if self.struct_type is None:
            self.struct_type = "name"
        if self.struct_type not in self._avail_structs:
            raise InputError(f"Invalid structure type: {self.struct_type}")
        
        # create vetee job object and read in coordinates
        job = get_struct_info(self.struct_input, self.struct_type)
        # if any of the following were given, compare to inputs
        # vetee obtained
        if self.charge is None:
            if job.charge is None:
                raise InputError("Unable to determine molecule charge.")
            self.charge = job.charge
        if self.multip is None:
            if job.multip is None:
                raise InputError("Unable to determine molecule multiplicity.")
            self.multip = job.multip  
        if job.charge != self.charge and job.charge is not None:
            print(f"Warning: default charge ({job.charge}) not equal \
                    \nto input charge ({self.charge}).")
        if job.multip != self.multip and job.multip is not None:
            print(f"Warning: default multiplicity ({job.multip}) not equal \
                    \nto input multiplicity ({self.multip}).")

        # check basis set and method agree if they were parsed in
        if self.struct_type == "com":
            basis = job._gaussian["gkeywords"]["basis"]
            method = job._gaussian["gkeywords"]["method"]
            if self.basis is None:
                assert basis is not None
                self.basis = basis
            elif self.basis != basis:
                print(f"Warning: parsed basis set ({basis}) not equal \
                        \nto input basis set ({self.basis}).")
            if self.method is None:
                assert method is not None
                self.method = method
            elif self.method != method:
                print(f"Warning: parsed method ({method}) not equal \
                        \nto input method ({self.method}).")

        # set default method if not read from com file
        if self.method is None:
            if self.prog == "openbabel":
                self.method = "mmff94"
            elif self.prog == "psi4":
                self.method = "hf"

        # write an xyz file to the current directory in order to generate
        # minimum dihedrals and openbabel object
        ofile = os.path.join(self.output_dir, "input_coords.xyz")
        write_xyz({"_coords": job.coords}, ofile)
        self.coords = job.xyz_coords
        
        # now determine minimum dihedrals by atom index
        self.min_diheds = get_min_dihed(ofile, self.charge, self.multip)
        self.obmol = create_obmol(ofile, self.charge, self.multip)
        
        # remove ring dihedrals where all four atoms are in a ring
        if self.no_ring_dihed:
            rings = get_rings(self.obmol)
            self.min_diheds = remove_ring_dihed(rings, self.min_diheds)
        self.num_dihed = len(self.min_diheds)
        # check that number of dihedral angles is not 0
        if self.num_dihed == 0:
            raise InputError("No dihedral angles could be generated for the input molecule.")

        # set atomic numbers
        self.atomic_nums = get_atomic_nums(self.obmol)


    def _set_output_dir(self):
        """Determine the name of the output directory.

        Parameters
        ----------
        structure : str
            Some identifier for the job. Example
            a string representing the name of the molecule.
        loc : str
            The parent directory to use as output.
            Defaults to "pwd", which means use the
            present working directory (current working
            directory). Another possible option is
            "home", which puts the output in the
            home directory (i.e. /user/home), but this
            option is only available for Linux users.
            If loc is not home or pwd, then the
            output will be generated in the given
            directory.

        Raises
        ------
        FileNotFoundError
            The user gave a location that does not exist.
        
        Notes
        -----
        The output is placed in kaplan_output under a job
        number, formatted as follows:
        loc/kaplan_output/job_0 # for the first job
        loc/kaplan_output/job_1 # for the second job
        etc.

        Returns
        -------
        output_dir : str
            The directory where the job output will
            be written.

        """
        if self.output_dir == "pwd":
            self.output_dir = os.path.join(os.getcwd(), "kaplan_output")
        elif self.output_dir == "home":
            self.output_dir = os.path.join(os.path.expanduser("~"), "kaplan_output")
        # don't make more subdirectories
        elif self.output_dir.endswith("kaplan_output"):
            pass
        else:
            self.output_dir = os.path.join(os.path.abspath(self.output_dir), "kaplan_output")

        # make an identifier for the molecule
        name = self.struct_input
        if self.struct_type in ("xyz", "com"):
            name = os.path.basename(name)
        # get rid of characters that should not be in file names
        # get rid of white space
        name = name.strip("!@#$%^&*()<>?/").replace(" ", "")

        # first check that there is a place to put the
        # output files
        if not os.path.isdir(self.output_dir):
            os.mkdir(self.output_dir)
            self.output_dir = os.path.join(self.output_dir, f"job_0_{name}")
            os.mkdir(self.output_dir)
            return
        # iterate over existing jobs to determine
        # dir_num for newest job
        dir_contents = os.scandir(self.output_dir)
        # keep track of how many jobs have been run
        dir_nums = []
        for val in dir_contents:
            val = val.name.split("_")
            if val[0] == "job" and len(val) == 3:
                try:
                    dir_nums.append(int(val[1]))
                except ValueError:
                    pass
        try:
            new_dir = f"job_{max(dir_nums)+1}_{name}"
        # max() arg is an empty sequence
        except ValueError:
            new_dir = f"job_0_{name}"
        self.output_dir = os.path.join(self.output_dir, new_dir)
        os.mkdir(self.output_dir)



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

        # only required input
        assert self.struct_input is not None

        # generate output directory
        self._set_output_dir()

        # read in structural information from file/openbabel/pubchem
        # writes an xyz file to the output directory
        self._update_geometry()

        # check that the inputs are valid
        self._check_input()


def read_input(input_file, new_output_dir=True):
    """Open a previously-written inputs pickle file.
    
    Parameters
    ----------
    input_file : str
        Explicit path and filename for inputs.pickle.
    new_output_dir : bool
        If True (default), then a new job directory will
        be generated for the inputs. If False, then the
        original directory (read from inputs object) will
        be used.

    Returns
    -------
    inputs object complete with obmol object and (possibly)
    update to output_dir.
    
    """
    inputs = Inputs()
    with open(input_file, "rb") as f:
        old_inputs = pickle.load(f)
    for var in old_inputs.__dict__:
        inputs.__dict__[var] = old_inputs.__dict__[var]
    # cannot pickle obmol object (since it's a swig object)
    # if the original input_coords.xyz file is not available,
    # have to regenerate obmol with the old coordinates
    # make an xyz file for openbabel to read
    orig_xyzfile = os.path.join(inputs.output_dir, "input_coords.xyz")

    if new_output_dir:
        # get rid of job_#_str directory
        # if the output_dir ends with a slash, then only the
        # slash is removed with dirname (could
        # nest kaplan_output by accident)
        if inputs.output_dir.endswith("/"):
            inputs.output_dir = inputs.output_dir[:-1]
        inputs.output_dir = os.path.dirname(inputs.output_dir)
        inputs._set_output_dir()

    if os.path.isfile(orig_xyzfile):
        obmol_obj = create_obmol(orig_xyzfile, inputs.charge, inputs.multip)
        if new_output_dir:
            os.rename(orig_xyzfile, os.path.join(inputs.output_dir, "input_coords.xyz"))

    else:
        full_coords = []
        for atom, coord in zip(inputs.atomic_nums, inputs.coords):
            full_coords.append([periodic_table(atom)] + list(coord))
        new_xyzfile = os.path.join(inputs.output_dir, "input_coords.xyz")
        write_xyz({"_coords": full_coords, "_comments": inputs.struct_input}, new_xyzfile)
        obmol_obj = create_obmol(new_xyzfile, inputs.charge, inputs.multip)

    # make sure the inputs are still valid
    inputs.obmol = obmol_obj
    inputs._check_input()
    return inputs
