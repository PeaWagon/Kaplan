
import numpy as np

# values for dihedral angles in degrees
MIN_VALUE = 0
MAX_VALUE = 360

class Pmem(object):
    """Population member of the ring."""

    def __init__(self, ring_loc, num_geoms, num_atoms,
                 current_mev, dihedrals=None):
        """Constructor for pmem object.

        Parameters
        ----------
        ring_loc : int
            Index of the ring where the pmem lives.
        num_geoms : int
            How many conformers we are trying to find.
        num_atoms : int
            Number of atoms in the input molecule. This
            determines the length of the dihedral angles
            list for each conformer.
        current_mev : int
            The mating event at which the pmem was constructed.
        dihedrals : list(list(int))
            The dihedrals for the pmem. Defaults to None.

        Attributes
        ----------
        dihedrals : np.array(shape(num_atoms-3, num_geoms),
                             dtype=int)
            List of integers representing the dihedral angles
            connecting the molecule under optimisation.

        Notes
        -----
        The birthday attribute keeps track of when the pmem
        was initialised. This will be helpful later in case
        the age of the pmems are pertinent to solving my
        problem.

        """
        self.ring_loc = ring_loc
        if dihedrals is None:
            # generate random dihedral angles (degrees)
            # each row is a set of dihedral angles for one conformer
            self.dihedrals = np.random.randint(MIN_VALUE, MAX_VALUE,
                             size=(num_geoms, num_atoms-3))
        else:
            self.dihedrals = dihedrals
        self.fitness = None
        self.birthday = current_mev

