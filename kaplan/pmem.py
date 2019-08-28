"""

This module contains the Pmem object.
Each instance of the pmem object keeps a record
of one set of conformers for a molecule.

A conformer is represented by a set of dihedral
angles. These dihedral angles are applied to
the initial geometry, which may or may not
be optimised (depending on the user's selection).

After the dihedral angles are applied, the fitness
is optimised using a forcefield to clean up
atom clashes and to explore the local energy space
for a better geometry. The pmem is scored based
on the energy for each conformer (each set of
dihedral angles) and the set of pairwise RMSD values.

"""

from math import factorial

import numpy as np

from kaplan.rmsd import calc_rmsd, apply_centroid
from kaplan.energy import run_energy_calc, MethodError, BasisSetError
# values for dihedral angles in radians
# convert radians to degrees for __str__ method
from kaplan.geometry import geometry_units,\
    update_obmol, set_coords, get_new_coordinates
from kaplan.inputs import Inputs, InputError
from kaplan.optimise import optimise_coords


class Pmem:
    """Population member of the ring."""

    def __init__(self, ring_loc, current_mev, num_geoms, num_diheds, dihedrals=None):
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
        num_diheds : int
            The number of dihedral angles to optimise.
        dihedrals : np.array(shape=(num_geoms, num_diheds),
                             dtype=float)
            The dihedrals for the pmem. Defaults to None.

        Raises
        ------
        AssertionError
            The dihedrals parameter does not have the correct
            shape (should be (num_geoms, num_diheds)).

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
        self.num_diheds = num_diheds
        if dihedrals is None:
            # generate random dihedral angles (degrees)
            # each row is a set of dihedral angles for one conformer
            inputs = Inputs()
            try:
                self.dihedrals = np.random.choice(
                    inputs.avail_diheds, replace=True,
                    size=(self.num_geoms, self.num_diheds)
                )
            # TypeError: 'NoneType' object cannot be interpreted as an integer
            except TypeError:
                raise InputError("Mising required inputs: num_geoms or num_diheds")
        else:
            assert dihedrals.shape == (self.num_geoms, self.num_diheds)
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

    def setup(self, major):
        """Setup a pmem for fitness evaluation.

        Parameters
        ----------
        major : bool
            If True, does a major optimisation for
            the pmem's conformers. If False, does a
            minor optimisation for the pmem's
            conformers.

        Notes
        -----
        This method is required if re-running a ring that
        contains pmems, where the inputs have been changed,
        such as program used for energy calculation.
        It is also required before a pmem can have its
        fitness evaluated.

        Coordinates are calculated as needed, but
        not stored. If the geometry specified by a
        set of dihedral angles does not converge,
        then the energy becomes None. For an rmsd
        calculation, if either of the geometries
        were non-convergent, the rmsd is None.
        This function should be called such that
        the fitness module can determine a pmem's
        fitness.

        Returns
        -------
        all_coords : list(np.array((num_atoms, 3), float))
            The length of the list is equal to num_geoms.
            This is a set of all optimised coordinates
            after the setup has complete.

        """
        all_coords = []
        for geom in range(self.num_geoms):
            all_coords.append(self._set_energy_get_coords(geom, major))

        if self.num_geoms == 1:
            return all_coords

        inputs = Inputs()

        for pair, (i, j) in enumerate(self.all_pairs_gen()):
            # make sure geometries are valid
            # need to explicitly check if None, because otherwise
            # numpy will complain that the check is ambiguous
            if all_coords[i] is not None and all_coords[j] is not None:
                if inputs.exclude_from_rmsd:
                    rmsd = calc_rmsd(
                        all_coords[i], all_coords[j],
                        inputs.atomic_nums, inputs.exclude_from_rmsd
                    )
                else:
                    rmsd = calc_rmsd(all_coords[i], all_coords[j])
                self.rmsds[pair] = [i, j, rmsd]
        return all_coords

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

    def _get_geometry(self, conf_index):
        """Get the base geometry after applying dihedral angles.

        Parameters
        ----------
        conf_index : int
            The index of the conformer for which to get
            the geometry.

        Notes
        -----
        There are two ways to get the geometry - using
        GOpt sets all of the dihedrals at once and tries
        to find a set of cartesian coordinates that
        correspond to these new internal coordinates.
        Openbabel applies a rotation matrix to all atoms
        on one side of each dihedral. Neither methods work
        for rings, and so constrained molecules are part
        of future work.

        When using Openbabel, this function has to set
        the obmol object to its original coordinates.
        This step is done to ensure all dihedrals are
        calculated relative to the same starting position.
        If a geometry cannot be applied using conformer
        n's dihedral angles, then None is returned.

        After getting these coordinates, the geometry
        should be locally optimised using a forcefield
        to account for clashing atoms. Then the geometry
        should be centred using the RMSD module.

        Returns
        -------
        coordinates : np.array((num_atoms, 3), float)
            The coordinates of the conformer after applying
            its dihedral angles. No centring is performed.
        None (if error occurs)

        """
        inputs = Inputs()

        # check if geometry optmisation is to be done using openbabel or GOpt
        if inputs.use_gopt:
            if "internal" not in inputs.extra:
                raise InputError(
                    "internal object required in inputs.extra dictionary."
                )
            # if GOpt doesn't converge, new_coords will be None
            new_coords = get_new_coordinates(
                inputs.extra["internal"], inputs.diheds, self.dihedrals[conf_index]
            )
            return new_coords

        if not hasattr(inputs, "obmol") or inputs.obmol is None:
            raise InputError("No input geometry has been set.")

        # reset obmol geometry to original coordinates
        set_coords(inputs.obmol, inputs.coords)
        try:
            result = update_obmol(inputs.obmol, inputs.diheds, self.dihedrals[conf_index])
        except AttributeError:
            raise InputError(f"Missing inputs. Unable to calculate geometry\
                               \nfor pmem {self.ring_loc}, conformer {conf_index}.")
        return result

    def _set_energy_get_coords(self, conf_index, major):
        """Set the energy for the input conformer index and get coordinates used.

        Parameters
        ----------
        conf_index : int
            The index of the conformer for which to get
            the coordinates and set the energy. Should
            be a value in the range 0 to num_geoms-1,
            inclusive.
        major : bool
            If False, do a minor optimisation. If True, the
            optimisation is major. Major optimisation is only
            performed for the initial geometry (if requested)
            and the final reporting of the best pmem at the end of
            evolution. Minor optimisation is performed during
            evolution, mostly to account for clashing atoms.

        Returns
        -------
        new_coords : np.array((num_atoms, 3), float)
            The coordinates after optimisation and centering,
            using the apply_centroid function from the RMSD
            module.
        None (if error occurs)

        """
        # get the coordinates after dihedral application
        raw_coords = self._get_geometry(conf_index)
        # if there was an issue, return None
        if raw_coords is None:
            return None

        # set energy of conformer
        try:
            self.energies[conf_index], new_coords = optimise_coords(raw_coords, major)
        # basis set and method errors should not be ignored,
        except (MethodError, BasisSetError, NotImplementedError, AssertionError):
            raise
        # if the exception is a psi4 exception, it is assumed
        # the energy calculation did not converge (for instance,
        # atoms being too close)
        # these types of errors should occur less often now that an
        # intermediate optimisation is being performed
        except Exception as e:
            # if there is a convergence error
            # give an energy of None
            # note: qcelemental.exceptions.ValidationError
            # is the atoms too close error
            print(e)
            print("Warning: psi4 did not converge.")
            self.energies[conf_index] = None
        new_coords = apply_centroid(new_coords)
        return new_coords
