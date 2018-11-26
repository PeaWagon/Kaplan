
from kaplan.pmem import Pmem
from kaplan.fitg import 
import numpy as np

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

class RingEmptyError(Exception):
    """Error that occurs when the ring has no pmems in it.

    Note
    ----
    This error will only be relevant once the ring
    extinction operators are enabled.

    """
    pass

class RingOverflowError(Exception):
    """Error that occurs when more pmems are added
    to the ring than there is space available."""
    pass

class Ring(object):
    """Data structure for genetic algorithm."""

    def __init__(self, num_geoms, num_atoms, num_slots,
                 pmem_dist, fit_form, coef_energy, coef_rmsd,
                 parser):
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

        Returns
        -------
        None

        """
        self.num_geoms = num_geoms
        self.num_atoms = num_atoms
        self.num_slots = num_slots
        self.pmem_dist = pmem_dist
        self.fit_form = fit_form
        if fit_form != 0:
            raise NotImplementedError("Only fit_form0 is available at this time.")
        self.coef_energy = coef_energy
        self.coef_rmsd = coef_rmsd
        self.parser = parser
        # make an empty ring
        self.num_filled = 0
        self.pmems = np.full(self.num_slots, None)

    def make_zmatrix(self, pmem_index):
        # generate a full geometry specification
        # in the form of a zmatrix using the dihedrals
        # found in the slot at pmem_index
        zmatrix = ""
        if self.pmems[pmem_index] == None:
            raise ValueError(f"Empty slot: {pmem_index}. No dihedrals with which to generate a zmatrix.")
        # generate zmatrix based on initial geom (parser)
        # replace the dihedrals with the dihedrals from the
        # pmem object
        # return a string
        return zmatrix

    def calc_fitness(self, pmem_index, zmatrix):
        """Calculate the fitness of a pmem.

        Parameters
        ----------
        pmem_index : int
            The location of the pmem in the ring.
        zmatrix : str
            The full geometry specification for
            the pmem at pmem_index.

        Notes
        -----
        Sets the value of pmem.fitness

        Returns
        -------
        None

        """
        pass

    def update(self, parent_index, child):
        """Add child to ring based on parent location.

        Parameters
        ----------
        parent_index : int
        children : pmem.dihedrals

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
        # check that adding pmems doesn't overflow ring
        num_avail = self.num_slots - self.num_filled
        try:
            assert num_avail >= num_pmems
        except AssertionError:
            raise RingOverflowError("Cannot add more pmems than space available in the ring.")
        # if self.num_filled == 0
        for i in range(self.num_filled, self.num_filled + num_pmems + 1):
            self.pmems[i] = Pmem(i, self.num_geoms,
                                 self.num_atoms, current_mev)
            # update pmem energies and fitness
        # if there are some pmems in the ring
        # they might not represent a contiguous segment
        # so go over each slot first and check that
        # it is empty
            

###############################################################

# attempt to enable ring indexing
# want to be able to do ring[30] and get the 30th pmem

#    def __iter__(self, i):
#        pass

#    def __getitem__(self, i):
#        for ind, pmem in enumerate(self.num_slots):
#            if ind > i:
#                
#        while i <

###############################################################    
    
