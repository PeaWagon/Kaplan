"""This module contains the Pmem object.
Each instance of the pmem object keeps a record
of one set of conformers for a molecule."""

import numpy as np

from saddle.errors import NotConvergeError

from kaplan.geometry import get_geom_from_dihedrals
from kaplan.fitg import calc_fitness

# values for dihedral angles in radians
MIN_VALUE = 0
MAX_VALUE = 2*np.pi


class Pmem:
    """Population member of the ring."""

    def __init__(self, ring_loc, current_mev, num_geoms, num_dihed, dihedrals=None):
        """Constructor for pmem object.

        Parameters
        ----------
        ring_loc : int
            Index of the ring where the pmem lives.
        current_mev : int
            The mating event at which the pmem was constructed.
        num_geoms : int
            How many conformers to search for.
        num_dihed : int
            The number of dihedral angles to optimise.
        dihedrals : np.array(shape=(num_geoms, num_dihed),
                             dtype=float)
            The dihedrals for the pmem. Defaults to None.

        Raises
        ------
        AssertionError
            The dihedrals parameter does not have the correct
            shape (should be (num_geoms, num_dihed)).

        Attributes
        ----------
        fitness : float
            The fitness of the pmem, as determined by
            the kaplan fitg module.
        birthday : int
            The mating event during which the pmem
            was generated.

        """
        self.ring_loc = ring_loc
        self.num_geoms = num_geoms
        self.num_dihed = num_dihed
        if dihedrals is None:
            # generate random dihedral angles (degrees)
            # each row is a set of dihedral angles for one conformer
            self.dihedrals = np.random.uniform(MIN_VALUE, MAX_VALUE,
                                               size=(self.num_geoms, self.num_dihed))
        else:
            assert dihedrals.shape == (self.num_geoms, self.num_dihed)
            self.dihedrals = dihedrals
        self.fitness = None
        self.birthday = current_mev


    def set_fitness(self):
        """Set the pmem's fitness.
        
        Notes
        -----
        Regarding non-convergence issues:
        Instead set the fitness to zero.
        Future work: determine fitness
        if 2 or more geometries converge
        (recalculate average and rmsds
        accordingly).
        
        """
        all_coords = []
        for geom in range(self.num_geoms):
            try:
                all_coords.append(get_geom_from_dihedrals(self.dihedrals[geom]))
            except NotConvergeError:
                self.fitness = 0
                return
        self.fitness = calc_fitness(all_coords)
