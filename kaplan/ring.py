"""Kaplan ring module.

This module is responsible for keeping track of
solution instances in the form of pmem objects.
This module is called by the tournament module
to perform updates to the population. This module
can also be used to determine summary statistics
about the population, such as mean, median,
and standard deviation for solution fitness,
energetic, and rmsd values. This module generates
an initial population at the start of evolution
by instantiating pmem objects.

The ring has a fixed number of slots, which are
iterable. Each slot either contains a pmem object
or None (when empty). The ring attribute num_filled
indicates how many slots are occupied at any given
time.

"""

from random import choice
from statistics import median, stdev, mean, StatisticsError

import numpy as np

from kaplan.inputs import Inputs, InputError
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
            return self.pmems[self.slot_num - 1]

    def __getitem__(self, key):
        """What happens when ring[integer] is called."""
        # need np.int64 since that is what goes into numpy arrays by default
        if not isinstance(key, int) and not isinstance(key, np.int64):
            raise KeyError("The ring cannot be indexed by non-integer values.")
        if key >= self.num_slots:
            raise KeyError("Given slot is larger than number of slots in ring.")
        return self.pmems[key]

    def __setitem__(self, key, value):
        """How to set ring[integer] = pmem.

        Notes
        -----
        Makes sure the pmem matches the current
        inputs when added to the ring, such as
        num_dihed, slot_num, and num_geoms.

        """
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
        # check pmem has correct inputs
        inputs = Inputs()
        assert inputs.num_geoms == value.num_geoms
        assert inputs.num_dihed == value.num_dihed
        # if not overwriting pmem slot, need to increment num_filled
        if self.pmems[key] is None:
            self.num_filled += 1
        self.pmems[key] = value

    def mating_radius(self, parent_index, mating_rad):
        """Return indices of the ring in the mating radius of parent.

        Parameters
        ----------
        parent_index : int
            Ring index for which to make a selection.
        mating_rad : int
            How many pmems to the left and right of parent
            are considered within the mating radius.

        Returns
        -------
        list(int)
            Indices representing ring slots that are
            in the mating radius of the parent slot.

        """
        # determine set of possible slots for the child to go
        # pick random index within +/-self.mating_rad of parent
        possible_slots = []
        # first check if the range loops round the ring
        if parent_index + mating_rad >= self.num_slots:
            possible_slots.extend(range(parent_index, self.num_slots))
            overflow = parent_index + mating_rad - self.num_slots
            possible_slots.extend(range(overflow + 1))
        else:
            possible_slots.extend(range(parent_index, parent_index + mating_rad + 1))
        # then check if loops around ring backwards
        # essentially checking if parent_index - mating_rad is negative
        if parent_index < mating_rad:
            possible_slots.extend(range(0, parent_index))
            backflow = self.num_slots - (mating_rad - parent_index)
            possible_slots.extend(range(backflow, self.num_slots))
        else:
            possible_slots.extend(range(parent_index - mating_rad, parent_index))

        # sometimes there are duplicates in the list due to
        # wrapping (i.e. mating_rad high and num_slots low)
        possible_slots = list(set(possible_slots))
        assert len(possible_slots) == 2 * mating_rad + 1

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
        # set the energetic and distance metrics needed for fitness calculation
        new_pmem.set_energy_rmsd()
        # note: the new pmem's energetic/rmsd values are not considered in the
        # normalisation. only pmems in the ring are considered
        # determine fitness value for the child
        # if normalise is False, then absolute fitness is used
        self.set_fitness(new_pmem)
        # check fitness vs current occupant (or empty slot)
        if self[potential_slot] is None or self[potential_slot].fitness is None or \
           self[potential_slot].fitness <= new_pmem.fitness:
            # add it there
            new_pmem._ring_loc = potential_slot
            self[potential_slot] = new_pmem

    def set_fitness(self, pmem):
        """Set the fitness for a given pmem.

        Parameters
        ----------
        pmem : kaplan.pmem.Pmem object
            The pmem whose fitness will be set.

        Notes
        -----
        This method was moved from the pmem module
        in order to accomodate for normalising
        energies and rmsd values. The normalisation
        is achieved by applying a z-score to each
        energy and rmsd. Normalising these values
        helps to make the coefficient inputs more
        intuitive, as energies and rmsd values are
        on a similar scale.

        z-score is defined as:
        z_score = (x-sample_mean)/sample_stdev

        The fitness of a pmem will change over the
        course of evolution (the sample mean and
        sample standard deviation will change), but
        its energy and rmsd values should remain
        constant. Therefore, the fitness of a pmem
        must be re-evaluated prior to:
        (1) tournament selection,
        (2) an extinction event (like plague
        or deluge) that orders pmems by fitness, and
        (3) the output-generation phase (when the
        best pmem is chosen).

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
        inputs = Inputs()
        if not inputs.normalise:
            pmem.set_fitness(inputs.fit_form, inputs.coef_energy, inputs.coef_rmsd)
            return pmem.fitness
        # make sure energies is not empty
        # rmsds can be empty in the case where there
        # is only one geometry per pmem (sum of empty list
        # is zero)
        # since we are maximising fitness, each energy value
        # should be negated (to favour low negative energies
        # and disfavour high positive energies - in the case
        # of forcefield calculations; quantum energies will
        # always be at least zero)
        valid_energies = [-e for e in pmem.energies if e is not None]
        valid_rmsds = [rmsd[2] for rmsd in pmem.rmsds if rmsd[2] is not None]
        # if all components for the fitness are None, return None
        if valid_energies == [] and valid_rmsds == []:
            pmem.fitness = None
            return None

        # account for statistics error in the case where there are 2 or
        # fewer entries for either energies or rmsds
        # mean requires 1 data point
        # stdev requires 2 data points
        try:
            mean_energy = -self.mean_energy
            stdev_energy = self.stdev_energy
        # division by zero
        except ZeroDivisionError:
            # there are no valid energies in the ring
            mean_energy = 0
            stdev_energy = 1
        # variance requires at least two data points
        except StatisticsError:
            stdev_energy = 1

        try:
            mean_rmsd = self.mean_rmsd
            stdev_rmsd = self.stdev_rmsd
        # division by zero
        except ZeroDivisionError:
            # there are no valid rmsds in the ring
            mean_rmsd = 0
            stdev_rmsd = 1
        # variance requires at least two data points
        except StatisticsError:
            stdev_rmsd = 1

        # normalise the energies and rmsds
        norm_energies = [(e - mean_energy) / stdev_energy for e in valid_energies]
        norm_rmsds = [(r - mean_rmsd) / stdev_rmsd for r in valid_rmsds]

        # now consider that the number of rmsd values and the number of energy
        # values are not the same
        energy_sum = sum(norm_energies) / inputs.num_geoms
        if inputs.num_geoms != 1:
            rmsd_sum = sum(norm_rmsds) / pmem.num_pairs
        else:
            rmsd_sum = 0

        if inputs.fit_form == 0:
            fitness = inputs.coef_energy * energy_sum + inputs.coef_rmsd * rmsd_sum
            pmem.fitness = fitness
            return fitness
        raise InputError("No such fitness formula.")

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

        # make sure ring has the same size as inputs
        assert inputs.num_slots == self.num_slots

        # check that adding pmems doesn't overflow ring
        num_avail = inputs.num_slots - self.num_filled
        try:
            assert num_avail >= num_pmems
        except AssertionError:
            raise RingOverflowError("Cannot add more pmems than space available in the ring.")

        # if there are no pmems in the ring, add a contiguous segment
        # and set fitness values
        # if setting normalised fitness, fill the slots first
        # then calculate fitness values
        if self.num_filled == 0:
            for i in range(0, num_pmems):
                self[i] = Pmem(i, current_mev, inputs.num_geoms,
                               inputs.num_dihed)
                self[i].set_energy_rmsd()
                if not inputs.normalise:
                    self[i].set_fitness(inputs.fit_form, inputs.coef_energy, inputs.coef_rmsd)
            if not inputs.normalise:
                return None
            for slot in range(0, num_pmems):
                self.set_fitness(self[slot])
            return None

        # if there are some pmems in the ring
        # they might not represent a contiguous segment
        # so go over each slot first and check that
        # it is empty
        total = self.num_filled + num_pmems
        for i in range(inputs.num_slots):
            if self.num_filled == total:
                if not inputs.normalise:
                    return None
                break
            if self[i] is None:
                self[i] = Pmem(i, current_mev, inputs.num_geoms,
                               inputs.num_dihed)
                self[i].set_energy_rmsd()
                if not inputs.normalise:
                    self[i].set_fitness(inputs.fit_form, inputs.coef_energy, inputs.coef_rmsd)
        for pmem in self.occupied:
            self.set_fitness(self[pmem])

    @property
    def occupied(self):
        """Return a list of slots that have a pmem."""
        slots = []
        for pmem in self:
            if pmem:
                slots.append(pmem.ring_loc)
        assert len(slots) == self.num_filled
        return slots

    @property
    def best_pmem(self):
        """Return the slot number and fitness value for the pmem with the highest fitness."""
        # need at least one pmem
        if not self.num_filled:
            raise RingEmptyError("No pmems in the ring.")
        best = None
        slot_num = None
        for slot in self.occupied:
            if self[slot].fitness is not None:
                if best is None or self[slot].fitness > best:
                    best = self[slot].fitness
                    slot_num = slot
        return (slot_num, best)

    @property
    def mean_fitness(self):
        return mean([self[i].fitness for i in self.occupied if self[i].fitness is not None])

    @property
    def mean_energy(self):
        total = 0
        count = 0
        for slot in self.occupied:
            for energy in self[slot].energies:
                # None means energy could not be calculated
                if energy is not None:
                    total += energy
                    count += 1
        return total / count

    @property
    def mean_rmsd(self):
        total = 0
        count = 0
        for slot in self.occupied:
            for rmsd in self[slot].rmsds:
                # None means one or two invalid geometries
                if rmsd[2] is not None:
                    total += rmsd[2]
                    count += 1
        return total / count

    @property
    def median_fitness(self):
        return median([self[i].fitness for i in self.occupied if self[i].fitness is not None])

    @property
    def median_energy(self):
        energy_vals = []
        for slot in self.occupied:
            for energy in self[slot].energies:
                # None means energy could not be calculated
                if energy is not None:
                    energy_vals.append(energy)
        return median(energy_vals)

    @property
    def median_rmsd(self):
        rmsd_vals = []
        for slot in self.occupied:
            for rmsd in self[slot].rmsds:
                # None means one or two invalid geometries
                if rmsd[2] is not None:
                    rmsd_vals.append(rmsd[2])
        return median(rmsd_vals)

    @property
    def stdev_fitness(self):
        return stdev([self[i].fitness for i in self.occupied if self[i].fitness is not None])

    @property
    def stdev_energy(self):
        energy_vals = []
        for slot in self.occupied:
            for energy in self[slot].energies:
                # None means energy could not be calculated
                if energy is not None:
                    energy_vals.append(energy)
        return stdev(energy_vals)

    @property
    def stdev_rmsd(self):
        rmsd_vals = []
        for slot in self.occupied:
            for rmsd in self[slot].rmsds:
                # None means one or two invalid geometries
                if rmsd[2] is not None:
                    rmsd_vals.append(rmsd[2])
        return stdev(rmsd_vals)
