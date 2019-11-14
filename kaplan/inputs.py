"""Default values for inputs as "borg"
class structure.

If the _options dictionary is changed
(by instantiating another of) Inputs
or DefaultInputs, the _options is updated
for all users.

This class is not responsible for checking
initial convergence or inputs related to the
energy calculations (i.e. it does not check
if the basis set and/or method are available
for the given program). These responsibilities
will be shifted to the energy module.

Parameters
----------
exclude_from_rmsd : list(int)
    Defaults to None, which means do not exclude
    any atomic numbers from the RMSD calculation.
    If not None, should be given as a list of
    integers, where each integer represents
    an atomic number to exclude from the RMSD
    calculation. For example:
    exclude = [1,8]
    would calculate the RMSD based on atoms
    that are not hydrogen or oxygen atoms.

"""


import os
import numpy as np
import pickle

from pybel import forcefields

from vetee.coordinates import CoordinatesError, pubchem_inchi_smiles, read_pubchem
from vetee.job import Job

from kaplan.geometry import update_obmol, create_obmol,\
    remove_ring_dihed, get_rings,\
    get_struct_info, get_atomic_nums, create_obmol_from_string,\
    get_coords, set_coords, write_coords,\
    get_torsions, filter_duplicate_diheds,\
    periodic_table

from kaplan.web import pubchem_request


# these values are used in energy calculations when no
# options are provided
# if a new program is added, add defaults for it
# how much RAM to use for psi4 calculations
# should be less than what your computer has available
hardware_inputs = {
    "psi4": {"RAM": "4 GB"},
    "openbabel": {},
    "other": {},
}


class InputError(Exception):
    """Raised when invalid input is given."""


class DefaultInputs:

    # programs available in kaplan
    _avail_progs = {"psi4", "openbabel"}

    # available structure inputs to use to generate coordinates
    # note: some require access to internet/requests library
    _avail_structs = {
        "sdf", "com", "xyz",        # file-based
        "smiles", "inchi",          # identifiers
        "name", "inchikey", "cid",  # requires pubchem
    }

    # available fitness function formats
    _avail_fit_form = [0]

    # characters that should not appear in name
    bad_chars = [" ", "~", "!", "@", "#", "$", "%", "^", "&", "*", "(", ")"]

    _options = {
        "output_dir": None,     # where to store the output
                                # if None, puts output in pwd
        "job_num": None,        # job number for output directory
        "name": None,           # used as a title for plots/directory names
                                # should not include any of the characters in bad_chars

        # mol inputs
        "prog": "openbabel",    # program to use to run energy calculations
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
        "min_dihed": True,      # select a minimum number of dihedral angles
        "exclude_from_rmsd": None,  # atomic numbers to exclude from RMSD calculations
                                    # if this value is not None, it should be list(int)
        "stop_at_conv": False,      # stop when the best pmem does not improve after
                                    # a certain number of mating events; False means run for
                                    # num_mevs, otherwise it should be an integer input specifying
                                    # how often to check the best pmem for improvement

        # the dihedral angles available that can be chosen for
        # new pmems, this should be a numpy array of values in
        # the range [-pi, pi)
        # by default this array is generated using numpy's
        # linspace function with 16 values
        # for example, if num is 4, then
        # the available dihedral angles (including those newly
        # generated and those mutated) would be:
        # -pi, -pi/2, 0, pi/2
        "avail_diheds": np.linspace(-np.pi, np.pi, num=16, endpoint=False),

        # geometry optimisation
        # uses the conjugate gradient method
        # major refers to initial and final geometries
        # minor refers to intermediate geometries (i.e. during evolution)
        # the only optimisation that can be turned off is the initial one
        "opt_init_geom": True,    # optimise the initial geometry using Openbabel forcefield
                                  # final/intermediate geometries will be optimised for all cases
        "major_tolerance": 1e-6,  # tolerance is the maximum difference in energy that is permitted
                                  # before the optimisation is considered to have converged
                                  # the tolerance may be set to a negative value, which means
                                  # that maxsteps will be completed (ignoring sampling)
        "major_maxsteps": 2500,   # if the tolerance condition is never satisfied, maxsteps is
                                  # how many iterations to complete in total
        "major_sampling": 100,    # sampling is how often to check the energy for the tolerance
                                  # in terms of number of conjugate gradient iterations
        "minor_tolerance": 0.1,
        "minor_maxsteps": 100,
        "minor_sampling": 10,

        # ea inputs
        "num_mevs": 100,
        "num_slots": 50,
        "init_popsize": 10,
        "mating_rad": 3,
        "num_geoms": 5,
        # default for num_swaps is equal to num_geoms//2 when not set
        "num_swaps": None,
        # default for num_muts is equal to (num_diheds*num_geoms)//5 when not set
        "num_muts": None,
        # default for num_cross is equal to num_geoms//2 when not set
        "num_cross": None,
        # crossover points is equal to num_diheds // 3 when not set
        # crossover points is set to 1 if num_diheds == 2
        # crossover is not completed if num_diheds == 1 (no cuts can
        # be made)
        # this parameter sets the maximum number of cuts to do
        # when doing a crossover mutation
        "max_cross_points": None,
        # where to do crossover; should be a list of values in the range
        # [1, num_diheds-1] inclusive
        # if None, crossover is done at any point in the list of dihedral angles
        "cross_points": None,

        # extinction operator inputs
        # each value is percent chance per mating event to apply operator
        "asteroid": 0.0,
        "plague": 0.0,
        "agathic": 0.0,
        "deluge": 0.0,

        # fitness function inputs
        "fit_form": 0,
        "coef_energy": 0.5,
        "coef_rmsd": 0.5,
        # normalise means that each energy and rmsd
        # value is normalised using z-score using all
        # available values in the ring
        # if normalise is false, the fitness is calculated
        # for each pmem as an absolute value
        "normalise": True,

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
        "num_diheds": None,
        # list of tuples, where each tuple has 4 integers
        # representing atom indices that make up a
        # dihedral angle
        "diheds": None,
        # extra is a dictionary that
        # contains keywords to use in the energy
        # calculations (such as convergence criteria)
        # these may be program specific
        "extra": {},
    }

    def __init__(self):
        self.__dict__ = self._options

    def __str__(self):
        return_str = ""
        for attr in self.__dict__:
            return_str += f"{attr}: {getattr(self, attr)}\n"
        return return_str


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
        self.output_dir = None
        self.job_num = None
        self.name = None
        self.prog = "openbabel"
        self.basis = "sto-3g"
        self.method = None
        self.struct_input = None
        self.struct_type = None
        self.charge = None
        self.multip = None
        self.no_ring_dihed = True
        self.min_dihed = True
        self.exclude_from_rmsd = None
        self.stop_at_conv = False
        self.avail_diheds = np.linspace(
            -np.pi, np.pi, num=16, endpoint=False
        )
        self.opt_init_geom = True
        self.major_tolerance = 1e-6
        self.major_maxsteps = 2500
        self.major_sampling = 100
        self.minor_tolerance = 0.1
        self.minor_maxsteps = 100
        self.minor_sampling = 10
        self.num_mevs = 100
        self.num_slots = 50
        self.init_popsize = 10
        self.mating_rad = 3
        self.num_geoms = 5
        self.num_swaps = None
        self.num_muts = None
        self.num_cross = None
        self.max_cross_points = None
        self.cross_points = None
        self.asteroid = 0.0
        self.plague = 0.0
        self.agathic = 0.0
        self.deluge = 0.0
        self.fit_form = 0
        self.coef_energy = 0.5
        self.coef_rmsd = 0.5
        self.normalise = True
        self.atomic_nums = None
        self.coords = None
        self.num_diheds = None
        self.diheds = None
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
            "max_cross_points",
            "num_diheds",
            "stop_at_conv",
            "major_maxsteps",
            "major_sampling",
            "minor_maxsteps",
            "minor_sampling",
        }
        expect_float = {
            "coef_energy", "coef_rmsd", "asteroid",
            "plague", "agathic", "deluge", "major_tolerance",
            "minor_tolerance",
        }
        expect_bool = {
            "no_ring_dihed", "normalise", "min_dihed",
            "opt_init_geom",
        }
        expect_list = {"exclude_from_rmsd", "cross_points"}
        expect_str = {"name"}

        # crossover does not make sense if there is only one dihedral angle
        if self.num_diheds == 1:
            self.num_cross = 0
            self.max_cross_points = 0

        # make sure all values have been set
        # and are set to the correct type
        for arg, val in self._options.items():
            if val is None:
                # defaults that are dependent on molecule size
                # or other inputs
                # val has to be reassigned otherwise
                # a ValueError occurs when int(None) is performed
                if arg == "num_swaps":
                    self.num_swaps = self.num_geoms // 2
                    val = self.num_swaps
                elif arg == "num_cross":
                    self.num_cross = self.num_geoms // 2
                    val = self.num_cross
                elif arg == "max_cross_points":
                    if self.num_diheds > 2:
                        self.max_cross_points = self.num_diheds // 3
                    else:
                        self.max_cross_points = 1
                    val = self.max_cross_points
                elif arg == "num_muts":
                    self.num_muts = (self.num_diheds * self.num_geoms) // 5
                    val = self.num_muts
                # these values are allowed to be None
                elif arg == "exclude_from_rmsd" or arg == "cross_points":
                    pass
                elif arg == "name":
                    self.name = self.get_name()
                    val = self.name
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
                except AssertionError:
                    raise InputError(
                        f"Value should be True or False and not a string: {arg} = {val}"
                    )
            elif arg in expect_list:
                try:
                    assert isinstance(val, list)
                except AssertionError:
                    if arg == "exclude_from_rmsd" and val is None:
                        pass
                    elif arg == "cross_points" and val is None:
                        pass
                    else:
                        raise InputError(
                            f"Value should be a list: {arg} = {val}"
                        )
            elif arg in expect_str:
                setattr(self, arg, str(val))

        assert isinstance(self.avail_diheds, np.ndarray)
        for value in self.avail_diheds:
            assert -np.pi <= value < np.pi
        assert 1 < len(self.avail_diheds)
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
        assert 5 < self.num_slots >= self.init_popsize
        assert self.num_mevs > 0
        assert self.num_geoms > 0
        assert 0 <= self.num_swaps <= self.num_geoms
        assert 0 <= self.num_cross <= self.num_geoms
        assert 0 <= self.num_muts <= self.num_diheds * self.num_geoms
        assert 2 <= self.mating_rad <= (self.num_slots - 1) // 2
        assert self.coef_energy >= 0
        assert self.coef_rmsd >= 0
        if self.exclude_from_rmsd is not None:
            for atom in self.exclude_from_rmsd:  # pylint: disable=not-an-iterable
                assert 1 <= atom <= 118
                assert isinstance(atom, int)
        assert 0 <= self.stop_at_conv
        assert 0 < self.num_diheds

        # check crossover-related items
        # crossover should not be called if there is only one dihedral
        if self.num_diheds == 1:
            assert self.num_cross == 0
        if 0 < self.num_cross:
            assert 1 <= self.max_cross_points <= self.num_diheds - 1
        if self.cross_points is not None:
            for value in self.cross_points:  # pylint: disable=not-an-iterable
                assert 1 <= value <= self.num_diheds - 1
                assert isinstance(value, int)

        # check the name for bad characters
        for char in self.bad_chars:
            try:
                assert char not in self.name
            except AssertionError:
                print("Do not put spaces in the name, use - or _ instead.")
                raise InputError(f"Bad character in name: '{char}'")
        try:
            if self.prog == "openbabel":
                assert self.method in forcefields
        except AssertionError:
            print("If running a quantum evaluation, make sure prog = \"psi4\".")
            raise InputError(f"Forcefield {self.method} not available in Openbabel.")
        assert 0 < self.major_maxsteps
        assert 0 < self.major_sampling
        assert 0 < self.minor_maxsteps
        assert 0 < self.minor_sampling

    def _update_geometry(self):
        """Update the OBMol object and the coordinates.

        Notes
        -----
        This method will now have two different approaches,
        one for smiles and one for everything else. Sometimes,
        Pubchem does not have the smiles string that is being
        input, but it is possible to create an OBMol object
        by reading the smiles string instead.

        This method adds the "pubchem_success" key to the
        self.extra dictionary. The value is True if the
        Pubchem query worked with the given input, and False
        otherwise. False also implies that Openbabel was used
        to generate the structure and not Pubchem.

        """
        # create vetee job object and read in coordinates
        try:
            job = get_struct_info(self.struct_input, self.struct_type)
            # check if molecule is flat
            # if flat, make it with Openbabel
            flat = all(np.allclose(atom[3], 0.0, atol=0.001) for atom in job.coords)
            if flat:
                print("Warning: flat molecule detected.")
                if self.struct_type not in ("com", "xyz"):
                    print("Changing input type to SMILES string.")
                    data = pubchem_inchi_smiles(self.struct_type, self.struct_input)
                    self.struct_input = data["_smiles_iso"]
                    self.struct_type = "smiles"
                    raise CoordinatesError
            self.extra["pubchem_success"] = True
        # vetee.coordinates.CoordinatesError occurs when
        # Pubchem database does not have what is given in the inputs
        except CoordinatesError as e:
            # can try to make OBMol from smiles/inchi, if given as input
            if self.struct_type not in ("smiles", "inchi"):
                raise e
            self.extra["pubchem_success"] = False
            print("Warning: Pubchem data ignored.")
            print("Generating OBMol and coordinates using Openbabel.")
            # make a new job object (from vetee)
            job = Job(self.struct_type, self.struct_input, "psi4")

            # try to read the smiles string with openbabel
            self.obmol = create_obmol_from_string(self.struct_type, self.struct_input)
            job.charge = self.obmol.GetTotalCharge()
            job.multip = self.obmol.GetTotalSpinMultiplicity()
            job._coords = []
            obmol_coords = get_coords(self.obmol)
            self.atomic_nums = get_atomic_nums(self.obmol)
            obmol_atoms = [periodic_table(a) for a in self.atomic_nums]
            for i, atom in enumerate(obmol_coords):
                job._coords.append([obmol_atoms[i], atom[0], atom[1], atom[2]])

        # if any of the following were given, compare to inputs
        # vetee obtained
        if self.charge is None:
            if job.charge is None:
                raise InputError("Unable to determine molecule charge.")
            self.charge = job.charge
        elif not isinstance(self.charge, int):
            raise InputError("Charge must be an integer.")

        if self.multip is None:
            if job.multip is None:
                raise InputError("Unable to determine molecule multiplicity.")
            self.multip = job.multip
        elif not isinstance(self.multip, int) or not self.multip >= 1:
            raise InputError("Multiplicity must be an integer value greater than 0.")

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
        write_coords(job.xyz_coords, job.atomic_nums, ofile)
        self.coords = job.xyz_coords

        # if pubchem succeeded in reading smiles /inchi string then
        # the obmol should already have been created
        if self.extra["pubchem_success"]:
            self.obmol = create_obmol(ofile, self.charge, self.multip)
            # set atomic numbers
            self.atomic_nums = get_atomic_nums(self.obmol)

        # now determine minimum dihedrals by atom index
        print("Selecting dihedral angles.")
        all_diheds = get_torsions(self.obmol)
        self.diheds = filter_duplicate_diheds(all_diheds, self.atomic_nums)

        # remove ring dihedrals where all four atoms are in a ring
        if self.no_ring_dihed:
            rings = get_rings(self.obmol)
            self.diheds = remove_ring_dihed(rings, self.diheds)
        self.num_diheds = len(self.diheds)
        # check that number of dihedral angles is not 0
        if self.num_diheds == 0:
            raise InputError("No dihedral angles could be generated for the input molecule.")

    def _set_output_dir(self):
        """Determine the path of the output directory.

        Notes
        -----
        The output is placed in kaplan_output under a job
        number, formatted as follows:
        ../kaplan_output/job_0 # for the first job
        ../kaplan_output/job_1 # for the second job
        etc.

        Returns
        -------
        output_dir : str
            The directory where the job output will
            be written.

        """
        if self.output_dir is None or self.output_dir == "":
            self.output_dir = os.path.join(os.getcwd(), "kaplan_output")
        else:
            # get rid of trailing slash
            if self.output_dir.endswith("/"):
                self.output_dir = self.output_dir[:-1]

            # don't make more subdirectories
            if os.path.basename(self.output_dir) == "kaplan_output":
                self.output_dir = os.path.abspath(self.output_dir)
            else:
                self.output_dir = os.path.join(os.path.abspath(self.output_dir), "kaplan_output")

            # make sure directory above kaplan_output exists
            if not os.path.isdir(os.path.dirname(self.output_dir)):
                raise InputError(f"No such directory exists: {self.output_dir}")

        # make an identifier for the molecule
        # should not be an issue if name is empty string
        if self.name is None:
            self.name = self.get_name()
        if any(c in self.name for c in self.bad_chars):
            print(f"Bad characters are: {self.bad_chars}")
            raise InputError(f"Name contains bad character: {self.name}")

        # first check that there is a place to put the
        # output files
        if not os.path.isdir(self.output_dir):
            os.mkdir(self.output_dir)

        # iterate over existing jobs to determine
        # dir_num for newest job
        dir_contents = os.scandir(self.output_dir)
        # keep track of how many jobs have been run
        dir_nums = []
        for val in dir_contents:
            val = val.name.split("_")
            if val[0] == "job" and len(val) >= 3:
                try:
                    dir_nums.append(int(val[1]))
                except ValueError:
                    pass
        try:
            job_num = max(dir_nums) + 1
        # max() arg is an empty sequence
        except ValueError:
            job_num = 0

        # zero-padding for job number means that ls will
        # list the jobs in numerical order (at least up until
        # 999999 jobs)
        new_dir = f"job_{job_num:06}_{self.name}"
        self.output_dir = os.path.join(self.output_dir, new_dir)
        os.mkdir(self.output_dir)

        return self.output_dir

    def write_inputs(self, fname="inputs.txt", overwrite=False):
        """Write the inputs object contents to a plain text file.

        Parameters
        ----------
        fname : str
            The name of the file to write the inputs to.
            Defaults to inputs.txt.
        overwrite : bool
            Defaults to False, which means an InputError
            will be raised if an attempt to overwrite
            an existing file is made. False means the file
            will be written regardless of whether an
            identical file name exists in the output
            directory.

        Raises
        ------
        FileExistsError
            The file already exists in the given directory.
        AssertionError
            The output_dir directory does not exist.

        Notes
        -----
        If the self.output_dir is not yet assigned,
        fname will be written to the current working
        directory.

        """
        assert os.path.isdir(self.output_dir)
        outname = os.path.join(self.output_dir, fname)
        try:
            # x means exclusive creation, failing if file exists
            with open(outname, "x") as f:
                f.write(self.__str__())
        except FileExistsError as e:
            if not overwrite:
                raise e
        with open(outname, "w") as f:
            f.write(self.__str__())

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
            # or a path designation
            if arg not in ("struct_input", "output_dir"):
                try:
                    val = val.lower()
                # value was integer or float
                except AttributeError:
                    pass
            # overwrite default value
            setattr(self, arg, val)

        # only required input
        assert self.struct_input is not None

        # make sure struct_type is valid otherwise
        # weird output directories are made
        if self.struct_type is None:
            self.struct_type = "name"
        elif self.struct_type not in self._avail_structs:
            raise InputError(f"Invalid structure type: {self.struct_type}")

        # generate output directory
        self._set_output_dir()

        # read in structural information from file/openbabel/pubchem
        # writes an xyz file to the output directory
        self._update_geometry()

        # check that the inputs are valid
        self._check_input()

        # write the contents to inputs.txt in the output_dir directory
        self.write_inputs()

        print("Done with inputs setup.")

    def get_name(self):
        """Returns a name for easy titles in plots and making a job directory."""

        if self.name is not None:
            return self.name

        # if the struct_input is an integer or float, it should be converted
        # to a string so that string-based methods work without error
        name = str(self.struct_input)

        # deal with files
        if self.struct_type in ("com", "xyz"):
            # strip the path and extension from the struct_input
            name = os.path.basename(os.path.splitext(name)[0])
            # get rid of white space
            # get rid of characters that should not be in file names
            for bad_char in self.bad_chars:
                name = name.replace(bad_char, "")

        elif self.struct_type == "name":
            name = name.replace(" ", "-")

        # get name for inchi or smiles since these will often have
        # special characters that cannot be used in directory names
        elif self.struct_type in ("smiles", "inchi"):
            try:
                data = read_pubchem(self.struct_type, self.struct_input)
                data = data["_synonyms"]
                if data == []:
                    raise CoordinatesError("No names available for input molecule.")
                name = data[0]
                for bad_char in self.bad_chars:
                    name = name.replace(bad_char, "")
            except CoordinatesError:
                name = "unknown_mol"

        return name

    def update_new(self, user_input):
        """New version of update_inputs using web module.

        Under construction**************

        """
        # can use openbabel
        if "struct_type" in ["com", "xyz", "sdf"]:
            pass
        elif "struct_type" in ["inchi", "smiles"]:
            pass
        # struct type is cid, name, inchikey
        # need to use pubchem
        else:
            pass

        # requests python library is needed to setup a
        # molecule structure with pubchem
        # if user_input["struct_type"] in ["name", "inchikey", "cid"]:
        #    try:
        #        pubchem_request(self.)
        #    if not REQUESTS:
        #        raise InputError("Requests library is required if the struct_type is name.")

        # pubchem_data = pubchem_request()


def read_input(job_inputs, new_output_dir=True):
    """Open a previously-written inputs pickle file or prep inputs object.

    Parameters
    ----------
    job_inputs : str or kaplan.inputs.Inputs object
        If a str, then it's an explicit path and
        filename for inputs.pickle. If object,
        it should have already been setup (except
        obmol).
    new_output_dir : bool, str
        If True (default), then a new job directory will
        be generated for the inputs. If False, then the
        original directory (read from inputs object) will
        be used. If new_output_dir is a str, then this will
        be the new directory for the inputs object (also
        the input_coords.xyz file will be written there).

    Notes
    -----
    If new_output_dir is not a boolean (True or False), then
    it should point to a Kaplan parent directory (example:
    /home/user/kaplan_output) or a directory that exists
    (example: /home/user). Do not put job_x_name in the
    output directory, because a new kaplan_output directory
    will then be added to said job_x_name directory (or an
    error will be raised in the event that said job directory
    does not exist).

    For jobs that are being restarted from a job that used
    a smiles string or inchi/inchikey as input, there is no
    guarantee that Openbabel generates the same atom
    ordering for the same string read in multiple times.
    Therefore it is recommended that users check that
    the structures are equivalent.

    Raises
    ------
    AssertionError
        The new_output_dir given does not exist.

    Returns
    -------
    inputs object complete with obmol object and (possibly)
    update to output_dir.

    """
    if isinstance(new_output_dir, str):
        assert os.path.isdir(new_output_dir)
    inputs = Inputs()

    if isinstance(job_inputs, str):
        with open(job_inputs, "rb") as f:
            old_inputs = pickle.load(f)
    else:
        # might have an issue if someone tries to pass DefaultInputs object
        assert isinstance(job_inputs, Inputs)
        old_inputs = job_inputs
    for var in old_inputs.__dict__:
        inputs.__dict__[var] = old_inputs.__dict__[var]
    # cannot pickle obmol object (since it's a swig object)
    # if the original input_coords.xyz file is not available,
    # have to regenerate obmol with the old coordinates
    # make an xyz file for openbabel to read
    orig_xyzfile = os.path.join(inputs.output_dir, "input_coords.xyz")

    # when set to True, still need old inputs dir as reference
    if new_output_dir is True:
        # get rid of job_#_str directory
        # if the output_dir ends with a slash, then only the
        # slash is removed with dirname (could
        # nest kaplan_output by accident)
        if inputs.output_dir.endswith("/"):
            inputs.output_dir = inputs.output_dir[:-1]
        # if the old output directory is being used, it will still have job_x_name
        inputs.output_dir = os.path.dirname(inputs.output_dir)
        inputs._set_output_dir()

    elif isinstance(new_output_dir, str):
        inputs.output_dir = new_output_dir
        inputs._set_output_dir()

    else:
        print("Warning: no new input_coords.xyz will be created.")
        print("The inputs object will still have the old directory as output_dir.")
        print(f"Any output for the current job will be written to:\n{inputs.output_dir}")

    # assume pubchem works if pubchem success key was never created
    try:
        pubchem_result = inputs.extra["pubchem_success"]
    except KeyError:
        pubchem_result = True

    if os.path.isfile(orig_xyzfile):
        if pubchem_result:
            obmol_obj = create_obmol(orig_xyzfile, inputs.charge, inputs.multip)
        else:
            obmol_obj = create_obmol_from_string(inputs.struct_type, inputs.struct_input)
            # make sure atomic numbers are in the same order
            # if this assertion fails, best to make a new inputs object
            # from scratch
            assert inputs.atomic_nums == get_atomic_nums(obmol_obj)
            # update the obmol_obj coordinates, since they are not the
            # same every time
            obmol_obj = set_coords(obmol_obj, inputs.coords)
        if new_output_dir:
            os.rename(orig_xyzfile, os.path.join(inputs.output_dir, "input_coords.xyz"))

    else:
        new_xyzfile = os.path.join(inputs.output_dir, "input_coords.xyz")
        write_coords(inputs.coords, inputs.atomic_nums, new_xyzfile)
        if pubchem_result:
            obmol_obj = create_obmol(new_xyzfile, inputs.charge, inputs.multip)
        else:
            obmol_obj = create_obmol_from_string(inputs.struct_type, inputs.struct_input)
            assert inputs.atomic_nums == get_atomic_nums(obmol_obj)
            obmol_obj = set_coords(obmol_obj, inputs.coords)

    # make sure the inputs are still valid
    inputs.obmol = obmol_obj
    inputs._check_input()
    return inputs


def get_latest_job(parent_dir):
    """Get the directory and job number for the last job in parent_dir.

    Parameters
    ----------
    parent_dir : str
        What is passed to inputs as output_dir. If this
        str does not end in kaplan_output, then said
        substring will be added.

    Returns
    -------
    tuple(str, int), where str is the directory name for the latest
    job and int is the last (corresponding) job number.

    """
    if not os.path.isdir(parent_dir):
        raise InputError(f"No such parent directory: {parent_dir}")

    if not parent_dir.endswith("kaplan_output"):
        parent_dir = os.path.join(parent_dir, "kaplan_output")

    # iterate over existing jobs to determine
    # dir_num for oldest job
    dir_contents = os.scandir(parent_dir)
    # keep track of how many jobs have been run
    max_job_num = -1
    max_job_name = None
    for val in dir_contents:
        val_list = val.name.split("_")
        if val_list[0] == "job" and len(val_list) == 3:
            try:
                job_num = int(val_list[1])
                if job_num > max_job_num:
                    max_job_num = job_num
                    max_job_name = val
            except ValueError:
                pass

    # kaplan_output exists but contains no job directories
    if max_job_name is None:
        raise InputError(f"No jobs in: {parent_dir}")

    # construct final output directory for last job
    return_dir = os.path.join(parent_dir, max_job_name)
    return (return_dir, max_job_num)
