"""
Consider the conversion of cartesian coordinates
to a z-matrix.

Likely to achieve using openbabel.

Need to clarify: does a z-matrix contain all
potential dihedral angles in a molecule?

Need to parse a coordinate system and generate
a list of dihedral angles.

Vetee already has the ability to read a z-matrix
file and write a z-matrix file using openbabel.
All that remains to be done is read a z-matrix
file and pull out dihedrals.

Generates an initial population to be used in evolution.

This will be an object file.

Important information that this object will contain:

* population size
* population structure (i.e. ring, koth, simple)
* molecule structure and connectivity
* instances of population members
*

"""

from random import choice

import numpy as np

from kaplan.inputs import Inputs
from kaplan.pmem import Pmem
from kaplan.fitg import sum_energies, sum_rmsds, calc_fitness
from kaplan.geometry import get_zmatrix_template, update_zmatrix, zmatrix_to_xyz


class RingEmptyError(Exception):
    """Error that occurs when the ring has no pmems in it.

    Note
    ----
    This error will only be relevant once the ring
    extinction operators are enabled.

    """


class RingOverflowError(Exception):
    """Error that occurs when more pmems are added
    to the ring than there is space available."""


class Ring:
    """Data structure for genetic algorithm."""

    def __init__(self):
        """Constructor for ring data structure.

        Parameters
        ----------
        num_geoms : int
            The number of geometries to optimize. Aka
            number of conformers to find for the input
            molecule. This number should not change.
        num_atoms : int
            The number of atoms that the input molecule
            has. This number should not change.
        num_slots : int
            The total number of slots in the ring that
            can be filled with population members (pmems).
            This number should not change.
        pmem_dist : int
            The distance that a pmem can be placed from
            the parent in number of slots.
        fit_form : int
            The number for the fitness function to use.
        coef_energy : float
            The coefficient for the sum of energies term
            for the fitness function.
        coef_rmsd : float
            The coefficient for the sum of rmsd values
            for the fitness function.
        parser : object
            The parser object from Vetee that contains
            information about molecular structure and
            energy calculations.

        Parameters
        ----------
        num_filled : int
            The population size contained within
            the ring. Should be less than or equal to the
            num_slots value. This variable is dynamic and
            changes depending on how many slots are currently
            filled in on the ring. Starts with a value of 0.
        pmems : np.ndarray(dtype=object)
            The list of pmem objects that are contained
            within the ring.
        zmatrix : str
            The original zmatrix specification (gzmat format)
            from the input geometry. Generated using geometry
            module (which uses openbabel).

        Returns
        -------
        None

        """
        inputs = Inputs()
        self.pmems = np.full(inputs.num_slots, None)
        # TODO: make sure zmatrix has charge and multip correctly set
        self.zmatrix = get_zmatrix_template(inputs.parser)

    def __getitem__(self, key):
        """What happens when ring[integer] is called."""
        if not isinstance(key, int):
            raise KeyError("The ring cannot be indexed by non-integer values.")
        inputs = Inputs()
        if key >= inputs.num_slots:
            raise KeyError("Given slot is larger than number of slots in ring.")
        return self.pmems[key]

    def __setitem__(self, key, value):
        """How to set ring[integer] = pmem."""
        if not isinstance(key, int):
            raise KeyError("The ring cannot be indexed by non-integer values.")
        inputs = Inputs()
        if key >= inputs.num_slots:
            raise KeyError("Given slot is larger than number of slots in ring.")
        if not isinstance(value, Pmem) and value is not None:
            raise KeyError("Ring should be filled with Pmem objects or None.")
        # in this case we are deleting a pmem
        if value is None:
            if self.pmems[key] is not None:
                inputs.num_filled -= 1
                self.pmems[key] = value
            return None
        # check that the pmem is being added to the same slot as ring_loc
        assert value.ring_loc == key
        # check that the pmem has the same num geoms and num atoms
        assert len(value.dihedrals) == inputs.num_geoms
        assert len(value.dihedrals[0]) == inputs.num_atoms - 3
        # if not overwriting pmem slot, need to increment num_filled
        if self.pmems[key] is None:
            inputs.num_filled += 1
        self.pmems[key] = value

    def set_fitness(self, pmem_index):
        """Set the fitness value for a pmem.

        Parameters
        ----------
        pmem_index : int
            The location of the pmem in the ring.

        Notes
        -----
        Sets the value of pmem.fitness

        Raises
        ------
        ValueError
            Slot is empty for given pmem_index.

        Returns
        -------
        None

        """
        if self.pmems[pmem_index] is None:
            raise ValueError(f"Empty slot: {pmem_index}.")
        inputs = Inputs()
        # construct zmatrices
        zmatrices = [update_zmatrix(self.zmatrix, self.pmems[pmem_index].dihedrals[i])
                     for i in range(inputs.num_geoms)]
        xyz_coords = [zmatrix_to_xyz(zmatrix) for zmatrix in zmatrices]
        # get fitness; set options for calculating energy
        energy = sum_energies(xyz_coords)
        rmsd = sum_rmsds(xyz_coords)
        fitness = calc_fitness(energy, rmsd)
        self.pmems[pmem_index].fitness = fitness

    def update(self, parent_index, child, current_mev):
        """Add child to ring based on parent location.

        Parameters
        ----------
        parent_index : int
            The location of the parent from which
            the child's location will be chosen.
        children : pmem.dihedrals
            Example, for num_geoms = 3, num_atoms = 8:
            [[3,4,5,3,4], [3,2,5,3,1], [6,3,2,5,6]]
        current_mev : int
            The current mating event. This method
            is called by the tournament. This value
            is used to set the birthday of new pmems
            (if the child is added to the ring).

        Notes
        -----
        This method is quite long. Perhaps it should
        be compartmentalised in the future.

        Returns
        -------
        None

        """
        inputs = Inputs()
        # determine fitness value for the child
        # construct zmatrices
        # TODO: since this is a repeat of the code from
        # set_fitness, it should be written as its own method
        zmatrices = [update_zmatrix(self.zmatrix, child[i]) for i in range(inputs.num_geoms)]
        xyz_coords = [zmatrix_to_xyz(zmatrix) for zmatrix in zmatrices]
        # get fitness
        energy = sum_energies(xyz_coords)
        rmsd = sum_rmsds(xyz_coords)
        fitness = calc_fitness(energy, rmsd)

        # TODO: see if this code should be replaced with negative
        # indices (since python lists are doubly-linked)
        # determine set of possible slots for the child to go
        # pick random index within +/-self.pmem_dist of parent
        possible_slots = []
        # first check if the range loops round the ring
        if parent_index + inputs.pmem_dist > inputs.num_slots:
            possible_slots.extend(range(parent_index, inputs.num_slots))
            overflow = parent_index + inputs.pmem_dist - inputs.num_slots
            possible_slots.extend(range(overflow+1))
        else:
            possible_slots.extend(range(parent_index, parent_index+inputs.pmem_dist+1))
        # then check if loops around ring backwards
        # essentially checking if parent_index - pmem_dist is negative
        if parent_index < inputs.pmem_dist:
            possible_slots.extend(range(0, parent_index))
            backflow = inputs.num_slots - (inputs.pmem_dist - parent_index)
            possible_slots.extend(range(backflow, inputs.num_slots))
        else:
            possible_slots.extend(range(parent_index-inputs.pmem_dist, parent_index))

        assert len(possible_slots) == 2*inputs.pmem_dist+1

        # select new child location
        chosen_slot = choice(possible_slots)
        # check fitness vs current occupant (or empty slot)
        if self[chosen_slot] is None or self[chosen_slot].fitness <= fitness:
            # add it there
            self[chosen_slot] = Pmem(chosen_slot, inputs.num_geoms,
                                     inputs.num_atoms, current_mev, child)
            self[chosen_slot].fitness = fitness

    def fill(self, num_pmems, current_mev):
        """Fill the ring with additional pmems.

        Notes
        -----
        This will fill a contiguous segment of the ring.
        Might have to change it to a try accept block
        where the pmems are added while there is space
        (otherwise I see this breaking for large t_size
        and small num_slots).

        Parameters
        ----------
        num_pmems : int
            Number of pmems to add to the ring.
        current_mev : int
            The current mating event (used to determine
            pmem age).

        Raises
        ------
        RingOverflowError
            Trying to add more pmems than slots available
            in ring.

        """
        inputs = Inputs()
        # check that adding pmems doesn't overflow ring
        num_avail = inputs.num_slots - inputs.num_filled
        try:
            assert num_avail >= num_pmems
        except AssertionError:
            raise RingOverflowError("Cannot add more pmems than space available in the ring.")
        # if there are no pmems in the ring, add a contiguous segment
        if inputs.num_filled == 0:
            for i in range(0, num_pmems):
                self.pmems[i] = Pmem(i, inputs.num_geoms,
                                     inputs.num_atoms, current_mev)
                self.set_fitness(i)
            inputs.num_filled += num_pmems
            return None
        # if there are some pmems in the ring
        # they might not represent a contiguous segment
        # so go over each slot first and check that
        # it is empty
        total = inputs.num_filled + num_pmems
        for i in range(inputs.num_slots):
            if inputs.num_filled == total:
                return None
            if self.pmems[i] is None:
                self.pmems[i] = Pmem(i, inputs.num_geoms, inputs.num_atoms,
                                     current_mev)
                self.set_fitness(i)
                inputs.num_filled += 1
