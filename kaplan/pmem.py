"""This module contains the Pmem object.
Each instance of the pmem object keeps a record
of one set of conformers for a molecule."""

from math import factorial

import numpy as np

from saddle.errors import NotConvergeError

from kaplan.geometry import get_geom_from_dihedrals
from kaplan.fitg import calc_fitness
from kaplan.rmsd import calc_rmsd
from kaplan.energy import run_energy_calc


# values for dihedral angles in radians
MIN_VALUE = 0
MAX_VALUE = 2*np.pi
# convert radians to degrees for __str__ method
RAD_TO_DEGREES = lambda x : x*180/np.pi


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
        all_coords : list
            Each member of the list is:
            np.array(shape=(n,3), dtype=float)
            Where n is the number of atoms. These
            arrays specify the xyz coordinates of
            one conformer. The length of the list
            is equal to num_geoms.
        energies : list
            List of energies for the conformers.
        rmsds: list
            Each member of the list is a tuple:
            conf_index1, conf_index2, rmsd
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
        self.all_coords = []
        self.energies = []
        self.rmsds = []
        self.fitness = None
        self.birthday = current_mev

    def __str__(self):
        """What happens when print(pmem) is called."""
        return_string = ""
        return_string += f"Slot #: {self.ring_loc}\n"
        for i, dihed_set in enumerate(self.dihedrals):
            return_string += f"\t Dihedrals {i}: "
            for dihed in dihed_set:
                return_string += f"{RAD_TO_DEGREES(float(dihed)):.2f} "
            return_string += "\n"
        return_string += f"\t Birthday: {self.birthday}\n"
        return_string += f"\t Fitness: {self.fitness}"
        return return_string

    def all_pairs_gen(self):
        """Yield indices of two geometries/conformers.

        Note
        ----
        This is a generator function.

        """
        for i in range(self.num_geoms-1):
            for j in range(i+1, self.num_geoms):
                yield (i, j)

    def set_fitness(self):
        """Set the pmem's fitness.

        Notes
        -----
        Also updates the all_coords attribute.
        If the dihedral angles do not converge
        to a geometry, then the geometry becomes
        None.
        Regarding non-convergence issues:
        If geometry is None, rmsd for that
        geometry (and any other geometry)
        is 0 and energy for that geometry
        is 0.
        
        """
        # get atomic coordinates for each conformer using dihedral angles
        for geom in range(self.num_geoms):
            try:
                self.all_coords.append(get_geom_from_dihedrals(self.dihedrals[geom]))
            except NotConvergeError:
                print("Warning: GOpt did not converge.")
                self.all_coords.append(None)

        # get energies of conformers
        for geometry in self.all_coords:
            if geometry is None:
                self.energies.append(0)
                continue
            try:
                self.energies.append(run_energy_calc(geometry))
            except Exception:
                # if there is a convergence error (atom too close)
                # give an energy of zero
                print("Warning: psi4 did not converge.")
                self.energies.append(0)

        # get rmsd for all pairs of conformers
        # n choose k = n!/(k!(n-k)!)
        num_pairs = int(factorial(self.num_geoms)/(2*factorial(self.num_geoms-2)))
        pairs = self.all_pairs_gen()
        for i in range(num_pairs):
            ind1, ind2 = next(pairs)
            if self.all_coords[ind1] is None or self.all_coords[ind2] is None:
                self.rmsds.append((ind1, ind2, 0))
            else:
                self.rmsds.append((ind1, ind2, calc_rmsd(self.all_coords[ind1],
                                                         self.all_coords[ind2])))
        
        # get fitness of pmem
        self.fitness = calc_fitness(self)
