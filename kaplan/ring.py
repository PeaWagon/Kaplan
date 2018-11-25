
from kaplan.pmem import Pmem

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


import numpy as np


class RingEmptyError(Exception):
    """Error that occurs when the ring has no pmems in it."""
    pass

class Ring(object):
    """Data structure for genetic algorithm."""

    def __init__(self, num_geoms, num_atoms, num_slots, 
                 num_filled, pmem_dist, mol_input_file):
#                 num_filled, dihedrals):
        """Constructor for ring data structure.

        Parameters
        ----------
        num_geoms : int
            The number of geometries to optimise. Aka
            number of conformers to find for the input
            molecule. This number should not change.
        num_atoms : int
            The number of atoms that the input molecule
            has. This number should not change.
        num_slots : int
            The total number of slots in the ring that
            can be filled with population members (pmems).
            This number should not change.
        num_filled : int
            The initial population size that goes into
            the ring. Should be less than or equal to the
            num_slots value. This variable is dynamic and
            changes depending on how many slots are currently
            filled in on the ring.
        pmem_dist : int
            The distance that a pmem can be placed from
            the parent in number of slots.

        Returns
        -------
        None

        """
        # my test
        # how many segments away from the parent a child can be
        # placed... will make this an input later
        self.max_distance = 2
        # make sure inputs are sound
        # TODO remove when complete
        try:
            assert num_slots >= num_filled
            assert isinstance(num_geoms, int)
            assert isinstance(num_atoms, int)
            assert isinstance(num_slots, int)
            assert isinstance(num_filled, int)

#            assert isinstance(dihedrals, list)
#            assert len(dihedrals) == num_atoms - 3
        except AssertionError as e:
            raise ValueError("Invalid input for Ring.")

#        self.dihedrals = dihedrals

        self.num_geoms = num_geoms
        self.num_atoms = num_atoms
        self.num_slots = num_slots
        self.num_filled = num_filled
        self.pmem_dist = pmem_dist
        self.mol_input_file = mol_input_file
        # Fill contiguous segment of ring with pmems.
        self.pmems = np.array([Pmem(i, num_geoms, mol_input_file) if i <= self.num_filled else None for i in range(self.num_slots)], dtype=object)
        # check that pmems have the same number of atoms as the mol_input_file
        assert self.num_atoms == self.pmems[0].num_atoms

    def update(self, parent_index, child):
        """Add child to ring based on parent location.

        Parameters
        ----------
        parent_index : int
        children : pmem

        Returns
        -------
        None

        """
        # NOTE: should check if person is added
        # if person is added, then increment self.num_filled


        # pick random index within +/-self.pmem_dist of parent
        # if slot is occupied: 
        #     compare fitness of child with current occupant
        #     if fitness no worse:
        #         replace occupant with child
        # else
        #     slot is empty so put child in slot
        pass

#    def update(self, child1, child2, loser1, loser2):
#        """All inputs are ring indices."""
#        # child1 and child2 have the same indices on
#        # the ring as their parents currently

#        # choose random location within the basket range
#        nest1 = np.random.choice(child1-self.max_distance, child1+self.max_distance)
#        nest2 = np.random.choice(child2-self.max_distance, child2+self.max_distance)
#        if ring.pmems[nest1] is not None:
#            del self.pmems[nest1]
#            del self.pmems[loser2]
#        self.pmems[loser1] = child1
#        self.pmems[loser2] = child2

    def fill(self, num_pmems):
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

        Raises
        ------
        AssertionError
            Trying to add more pmems than slots available
            in ring.

        """
        # check that adding pmems doesn't overflow ring
        num_avail = self.num_slots - self.num_filled
        assert num_avail >= num_pmems
        for i in range(self.num_filled, self.num_filled + num_pmems + 1):
            self.pmems[i] = Pmem(i, self.num_geoms, self.mol_input_file)
            

###############################################################

# attempt to enable ring indexing

#    def __iter__(self, i):
#        pass

#    def __getitem__(self, i):
#        for ind, pmem in enumerate(self.num_slots):
#            if ind > i:
#                
#        while i <





###############################################################
    def count_occupied(self):
        """Determine number of filled slots in the ring."""
        filled = 0
        for pmem in self.pmems:
            if pmem is not None:
                filled += 1
        return filled


    
    
