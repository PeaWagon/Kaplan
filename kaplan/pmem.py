"""This module contains the Pmem object.
Each instance of the pmem object keeps a record
of one set of conformers for a molecule."""

from math import factorial

import numpy as np

from kaplan.rmsd import calc_rmsd, apply_centroid
from kaplan.energy import run_energy_calc, MethodError, BasisSetError
# values for dihedral angles in radians
# convert radians to degrees for __str__ method
from kaplan.geometry import MIN_VALUE, MAX_VALUE, geometry_units, update_obmol, set_coords
from kaplan.inputs import Inputs, InputError


class Pmem:
    """Population member of the ring."""

    def __init__(self, ring_loc, current_mev, num_geoms, num_dihed, dihedrals=None):
        """Constructor for pmem object.

        Parameters
        ----------
        ring_loc : int, None
            Index of the ring where the pmem lives.
            Is None if the pmem is generated but not
            yet in the ring.
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
            Starts as list of None. If an energy can't be
            calculated, then it retains a value of None.
        rmsds: list
            Each member of the list is a list containing:
            conf_index1, conf_index2, rmsd. The rmsd
            starts as None. If a geometry cannot
            be made for conformer 1 and/or conformer
            2, then the rmsd is kept as None. Rmsd
            can still be calculated even if the geometry
            does not converge; the calculation only
            requires a geometry to be made.
        fitness : float
            The fitness of the pmem, as set by the ring
            module.
        birthday : int
            The mating event during which the pmem
            was generated.

        """
        self._ring_loc = ring_loc
        self.num_geoms = num_geoms
        self.num_dihed = num_dihed
        if dihedrals is None:
            # generate random dihedral angles (degrees)
            # each row is a set of dihedral angles for one conformer
            try:
                self.dihedrals = np.random.uniform(
                    MIN_VALUE, MAX_VALUE, size=(self.num_geoms, self.num_dihed)
                )
            # TypeError: 'NoneType' object cannot be interpreted as an integer
            except TypeError:
                raise InputError("Mising required inputs: num_geoms or num_dihed")
        else:
            assert dihedrals.shape == (self.num_geoms, self.num_dihed)
            self.dihedrals = dihedrals
        self.energies = [None for _ in range(self.num_geoms)]
        self.rmsds = [[i, j, None] for i, j in self.all_pairs_gen()]
        self.fitness = None
        self._birthday = current_mev

    def __str__(self):
        """What happens when print(pmem) is called.

        Notes
        -----
        Dihedral angles are printed in degrees.

        """
        return_string = ""
        return_string += f"Slot #: {self.ring_loc}\n"
        for i, dihed_set in enumerate(self.dihedrals):
            return_string += f"\t Dihedrals {i}: "
            for dihed in dihed_set:
                return_string += f"{geometry_units['radians']['degrees'](float(dihed)):.2f} "
            return_string += "\n"
        return_string += f"\t Birthday: {self.birthday}\n"
        return_string += f"\t Fitness: {self.fitness}\n"
        return_string += f"\t Energies: {self.energies}\n"
        return_string += f"\t RMSDs: {self.rmsds}"
        return return_string

    def __iter__(self):
        """Called at the start of iteration."""
        self.conf_num = 0
        return self

    def __next__(self):
        """Called for iteration over conformers.

        Notes
        -----
        Dihedrals give here are in radians.

        """
        if self.conf_num == self.num_geoms:
            raise StopIteration
        else:
            self.conf_num += 1
            return self.dihedrals[self.conf_num - 1]

    def all_pairs_gen(self):
        """Yield indices of two geometries/conformers.

        Note
        ----
        This is a generator function.

        """
        for i in range(self.num_geoms - 1):
            for j in range(i + 1, self.num_geoms):
                yield (i, j)

    @property
    def num_pairs(self):
        """How many ways we can pick two conformers (ignoring order)."""
        if self.num_geoms == 1:
            return 0
        # n choose k = n!/(k!(n-k)!)
        return int(factorial(self.num_geoms) / (2 * factorial(self.num_geoms - 2)))

    @property
    def ring_loc(self):
        return self._ring_loc

    @property
    def birthday(self):
        return self._birthday

    def get_geometry(self, conf_index):
        """Get a centred geometry for a given conformer index.

        Parameters
        ----------
        conf_index : int
            The index of the conformer for which to get
            the geometry.

        Notes
        -----
        First this function has to set the obmol object
        to its original coordinates. This step is done
        to ensure all dihedrals are calculated relative
        to the same starting position. If a geometry
        cannot be applied using conformer n's dihedral
        angles, then None is returned.

        Returns
        -------
        coordinates : np.array((num_atoms, 3), float)
            The coordinates of the conformer, after centering
            using apply_centroid function from rmsd module.
        None (if error occurs)

        """
        inputs = Inputs()
        if not hasattr(inputs, "obmol") or inputs.obmol is None:
            raise InputError("No input geometry has been set.")
        # reset obmol geometry to original coordinates
        set_coords(inputs.obmol, inputs.coords)
        try:
            result = update_obmol(inputs.obmol, inputs.min_diheds, self.dihedrals[conf_index])
        except AttributeError:
            raise InputError(f"Missing inputs. Unable to calculate geometry\
                               \nfor pmem {self.ring_loc}, conformer {conf_index}.")
        return apply_centroid(result)

    def set_energy_get_coords(self, conf_index):
        """Set the energy for the input conformer index and get coordinates used."""
        new_coords = self.get_geometry(conf_index)
        # if there was an issue, return None
        if new_coords is None:
            return None
        # set energy of conformer
        try:
            self.energies[conf_index] = run_energy_calc(new_coords)
        # basis set and method errors should not be ignored,
        except (MethodError, BasisSetError):
            raise
        # if the exception is a psi4 exception, it is assumed
        # the energy calculation did not converge (for instance,
        # atoms being too close)
        except Exception as e:
            # if there is a convergence error
            # give an energy of None
            # note: qcelemental.exceptions.ValidationError
            # is the atoms too close error
            print(e)
            print("Warning: psi4 did not converge.")
            self.energies[conf_index] = None
        return new_coords

    def set_energy_rmsd(self):
        """Set the pmem's fitness precursors, energies and rmsds.

        Notes
        -----
        Coordinates are calculated as needed, but
        not stored. If the geometry specified by a
        set of dihedral angles does not converge,
        then the energy becomes None. For an rmsd
        calculation, if either of the geometries
        were non-convergent, the rmsd is None.
        This function should be called such that
        the ring can determine a pmem's fitness.

        """
        # make sure everything is empty before setting up fitness
        assert self.energies == [None for _ in range(self.num_geoms)]
        assert self.rmsds == [[i, j, None] for i, j in self.all_pairs_gen()]

        # if there is only one geometry, fitness should only be based on
        # the energy (since no RMSD can be calculated)
        if self.num_geoms == 1:
            self.set_energy_get_coords(0)
            return None

        # keep track of which geometries we have already calculated centroid for
        # 0 different from None is needed in case geometry fails
        coords = [0 for _ in range(self.num_geoms)]  # retain coords for this function only

        # get rmsd for all pairs of conformers and energies for all geoms
        for i, (ind1, ind2, _) in enumerate(self.rmsds):
            if coords[ind1] is 0:
                coords[ind1] = self.set_energy_get_coords(ind1)
            if coords[ind2] is 0:
                coords[ind2] = self.set_energy_get_coords(ind2)

            # if either of the geometries could not be constructed,
            # keep the rmsd as None
            if coords[ind1] is not None and coords[ind2] is not None:
                self.rmsds[i][2] = calc_rmsd(coords[ind1], coords[ind2])

    def set_fitness(self, fit_form, coef_energy, coef_rmsd):
        """Calculate the absolute fitness of a pmem.

        Parameters
        ----------
        fit_form : int
            Fitness formula to use (from inputs).
        coef_energy : float
            Energy term coefficient (from inputs).
        coef_rmsd : float
            Root-mean-square deviation coefficient
            (from inputs).

        Notes
        -----
        This method calculates the fitness of a
        population member (pmem). The fitness is
        broken into two main parts: energy and
        rmsd. The energy is calculated using a
        user-specified quantum chemical method
        and basis set for each geometry in the pmem.
        The rmsd is calculated as all the possible
        pairs of rmsd between geometries.

        fit_form : int
            Represents the fitness formula to use.
            The only value currently available is 0,
            where fitness = CE*SE + Crmsd*Srmsd.

        Returns
        -------
        fitness : float

        """
        # make sure energies is not empty
        # rmsds can be empty in the case where there
        # is only one geometry per pmem (sum of empty list
        # is zero)
        valid_energies = [e for e in self.energies if e is not None]
        valid_rmsds = [rmsd[2] for rmsd in self.rmsds if rmsd[2] is not None]
        # if all components for the fitness are None, return None
        if valid_energies == [] and valid_rmsds == []:
            self.fitness = None
            return None
        if fit_form == 0:
            fitness = -1 * coef_energy * sum(valid_energies) + coef_rmsd * sum(valid_rmsds)
            self.fitness = fitness
            return fitness
        raise InputError("No such fitness formula.")
