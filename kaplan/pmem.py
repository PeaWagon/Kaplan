
from kaplan.fitg import get_fitness
from kaplan.geometry import get_parser

class Pmem(object):
    """Population member of the ring."""

    def __init__(self, ring_loc, num_geoms, mol_input_file):
        """Constructor for pmem object.

        Parameters
        ----------
        ring_loc : int
            Index of the ring where the pmem lives.
        num_geoms : int
            How many conformers we are trying to find.
        mol_input_file : str
            The name of the file representing the
            geometry specifications.

        Attributes
        ----------
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
