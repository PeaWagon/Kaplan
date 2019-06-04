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

    def __init__(self, num_slots, init_popsize):
        """Constructor for ring data structure.

        Parameters
        ----------
        num_slots : int
            The total number of slots in the ring that
            can be filled with population members (pmems).
            This number should not change.
        init_popsize : int
            The number of pmems for the initial population.

        Attributes
        ----------
        num_filled : int
            The population size contained within
            the ring. Should be less than or equal to the
            num_slots value. This variable is dynamic and
            changes depending on how many slots are currently
            filled in on the ring. Starts with a value equal
            to the init_popsize parameter (from inputs).
        pmems : np.array(dtype=object)
            The array of pmem objects that are contained
            within the ring.

        Returns
        -------
        None

        """
        self.num_slots = num_slots
        self.pmems = np.full(self.num_slots, None)
        # create initial population
        self.num_filled = 0
        self.fill(init_popsize, 0)

    def __str__(self):
        """What happens when print(ring) is called."""
        return_string = f"Total ring slots: {self.num_slots}\n"
        return_string += f"Filled ring slots: {self.num_filled}\n"
        return_string += f"Available ring slots: {self.num_slots-self.num_filled}\n"
        for i, pmem in enumerate(self):
            if pmem:
                return_string += str(pmem) + "\n"
            else:
                return_string += f"Slot #: {i} => Empty\n"
        return return_string

    def __iter__(self):
        """Called at the start of iteration."""
        self.slot_num = 0
        return self
    
    def __next__(self):
        """Called for iteration over slots."""
        if self.slot_num == self.num_slots:
            raise StopIteration
        else:
            self.slot_num += 1
            return self.pmems[self.slot_num-1]

    def __getitem__(self, key):
        """What happens when ring[integer] is called."""
        # need np.int64 since that is what goes into numpy arrays by default
        if not isinstance(key, int) and not isinstance(key, np.int64):
            raise KeyError("The ring cannot be indexed by non-integer values.")
        if key >= self.num_slots:
            raise KeyError("Given slot is larger than number of slots in ring.")
        return self.pmems[key]

    def __setitem__(self, key, value):
        """How to set ring[integer] = pmem."""
        if not isinstance(key, int):
            raise KeyError("The ring cannot be indexed by non-integer values.")
        if key >= self.num_slots:
            raise KeyError("Given slot is larger than number of slots in ring.")
        if not isinstance(value, Pmem) and value is not None:
            raise KeyError("Ring should be filled with Pmem objects or None.")
        # in this case we are deleting a pmem
        if value is None:
            if self.pmems[key] is not None:
                self.num_filled -= 1
                self.pmems[key] = value
            return None
        # check that the pmem is being added to the same slot as ring_loc
        assert value.ring_loc == key
        # if not overwriting pmem slot, need to increment num_filled
        if self.pmems[key] is None:
            self.num_filled += 1
        self.pmems[key] = value


    def mating_radius(self, parent_index, pmem_dist):
        """Return indices of the ring in the mating radius of parent.

        Parameters
        ----------
        parent_index : int
            Ring index for which to make a selection.
        pmem_dist : int
            How many pmems to the left and right of parent
            are considered within the mating radius.
        
        Returns
        -------
        list(int)
            Indices representing ring slots that are
            in the mating radius of the parent slot.
        
        """
        # TODO: see if this code should be replaced with negative
        # indices (since python lists are doubly-linked)
        # determine set of possible slots for the child to go
        # pick random index within +/-self.pmem_dist of parent
        possible_slots = []
        # first check if the range loops round the ring
        if parent_index + pmem_dist >= self.num_slots:
            possible_slots.extend(range(parent_index, self.num_slots))
            overflow = parent_index + pmem_dist - self.num_slots
            possible_slots.extend(range(overflow+1))
        else:
            possible_slots.extend(range(parent_index, parent_index+pmem_dist+1))
        # then check if loops around ring backwards
        # essentially checking if parent_index - pmem_dist is negative
        if parent_index < pmem_dist:
            possible_slots.extend(range(0, parent_index))
            backflow = self.num_slots - (pmem_dist - parent_index)
            possible_slots.extend(range(backflow, self.num_slots))
        else:
            possible_slots.extend(range(parent_index-pmem_dist, parent_index))

        # sometimes there are duplicates in the list due to
        # wrapping (i.e. pmem_dist high and num_slots low)
        possible_slots = list(set(possible_slots))
        assert len(possible_slots) == 2*pmem_dist+1

        return possible_slots


    def update(self, child, potential_slot, current_mev):
        """Add child to ring based on parent location.

        Parameters
        ----------
        child : np.array(shape=(num_geoms, num_dihed),
                         dtype=float)
        potential_slot : int
            The slot to place the child.
        current_mev : int
            The current mating event. This method
            is called by the tournament. This value
            is used to set the birthday of new pmems
            (if the child is added to the ring).

        Returns
        -------
        None

        """
        num_geoms, num_dihed = child.shape
        # create new pmem object using dihedral angles
        new_pmem = Pmem(None, current_mev, num_geoms, num_dihed, dihedrals=child)
        # determine fitness value for the child
        new_pmem.set_fitness()
        # check fitness vs current occupant (or empty slot)
        if self[potential_slot] is None or self[potential_slot].fitness <= new_pmem.fitness:
            # add it there
            new_pmem.ring_loc = potential_slot
            self[potential_slot] = new_pmem


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
        # generate necessary inputs
        inputs = Inputs()
        
        # check that adding pmems doesn't overflow ring
        num_avail = inputs.num_slots - self.num_filled
        try:
            assert num_avail >= num_pmems
        except AssertionError:
            raise RingOverflowError("Cannot add more pmems than space available in the ring.")
        
        # if there are no pmems in the ring, add a contiguous segment
        # and set fitness values
        if self.num_filled == 0:
            for i in range(0, num_pmems):
                self[i] = Pmem(i, current_mev, inputs.num_geoms,
                               inputs.num_dihed)
                self[i].set_fitness()
            return None
        
        # if there are some pmems in the ring
        # they might not represent a contiguous segment
        # so go over each slot first and check that
        # it is empty
        total = self.num_filled + num_pmems
        for i in range(inputs.num_slots):
            if self.num_filled == total:
                return None
            if self[i] is None:
                self[i] = Pmem(i, current_mev, inputs.num_geoms,
                               inputs.num_dihed)
                self[i].set_fitness()
