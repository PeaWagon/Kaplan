

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



def get_dihedrals(input_file):
    with open(input_file, 'r') as f:
        pass
        
    
    




class Ring(object):
    """Data structure for genetic algorithm."""

    def __init__(self, num_geoms, num_atoms, num_slots, 
                 num_filled):
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
        # Fill contiguous segment of ring with pmems.
        self.pmems = np.array([Pmem(i, self.num_geoms, self.num_atoms) if i <= self.num_filled else None for i in range(self.num_slots)], dtype=object)

    def update(self, child1, child2, loser1, loser2):
        """All inputs are ring indices."""
        # child1 and child2 have the same indices on
        # the ring as their parents currently

        # choose random location within the basket range
        nest1 = np.random.choice(child1-self.max_distance, child1+self.max_distance)
        nest2 = np.random.choice(child2-self.max_distance, child2+self.max_distance)
        if ring.pmems[nest1] is not None:
            del self.pmems[nest1]
            del self.pmems[loser2]
        self.pmems[loser1] = child1
        self.pmems[loser2] = child2
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

class Pmem(Ring):
    """Population member of the ring."""

    def __init__(self, ring_loc, num_geoms, num_atoms):
        """

        Parameters
        ----------
        ring_loc : int
            Index of the ring where the pmem lives.
        dihedrals : np.array(shape(num_atoms-3, num_geoms),
                             dtype=int)
            List of integers representing the dihedral angles
            connecting the molecule under optimisation.

        """
        self.ring_loc = ring_loc
        self.num_geoms = num_geoms
        self.num_atoms = num_atoms
#        self._fitness = None
        # generate random dihedral angles (degrees)
        # each row is a set of dihedral angles for one conformer
        self.dihedrals = np.random.randint(0, 359,
                         size=(self.num_geoms,self.num_atoms-3))
        # TODO change to actually calculate the energy
        self.energies = np.zeros(shape=(self.num_geoms,),
                                 dtype=float)
#        self._fitness = np.zeros(shape=(self.num_slots,), dtype=float)


#    def is_occupant(self):
#        """Determine if Pmem is part of Ring or just a placeholder."""
#        if self.ring_loc is not None:
#            return True
#        else:
#            return False

#    @fitness.setter
#    def fitness(self, value):
#        assert isinstance(value, float)
#        assert value >= 0
#        self._fitness = value
    
    
