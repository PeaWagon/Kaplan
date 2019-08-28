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

The mean, median, and standard deviation can be
calculated for fitness, energy, and RMSD. If these
values cannot be calculated, the default behaviour is
to accept the StatisticsError (from statistics module)
and return None. The main case where these cannot
be calculated is the event where only one pmem
is in the ring. Also errors can occur if energies
cannot be calculated for all geometries in the ring.

You need at least 2 data points for stdev calculation.
The StatisticsError can occur if the output module
is called right after an extinction event (i.e. only
one pmem left).

This module is not responsible for updating or setting
the fitness values of any pmems.

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
        num_diheds, slot_num, and num_geoms.

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
        assert inputs.num_diheds == value.num_diheds
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

    def update(self, child, potential_slot):
        """Add child to ring based on potential_slot location.

        Parameters
        ----------
        child : kaplan.pmem.Pmem object
            A pmem that wants to be in the ring. Should
            have its fitness already set.
        potential_slot : int
            The slot to place the child. If this slot
            is occupied, its fitness should be set.

        Returns
        -------
        success : bool
            True if new pmem was added to the ring. False
            if new pmem was not added to the ring.

        """
        # in the rare event that the child has no valid
        # energies or RMSD values, its fitness will be None
        if child.fitness is None:
            print("Warning. Child fitness was None.")
            if self[potential_slot] is None:
                success = True
            else:
                success = False

        # check fitness vs current occupant (or empty slot)
        elif self[potential_slot] is None or \
                self[potential_slot].fitness is None or \
                self[potential_slot].fitness <= child.fitness:
            success = True

        # fitness of new pmem is not as good as old pmem
        else:
            success = False

        if success:
            child._ring_loc = potential_slot
            self[potential_slot] = child
        return success

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
        # don't calculate the fitness for any pmems
        if self.num_filled == 0:
            for i in range(num_pmems):
                self[i] = Pmem(i, current_mev, inputs.num_geoms,
                               inputs.num_diheds)
                self[i].setup(major=False)
            return None

        # if there are some pmems in the ring
        # they might not represent a contiguous segment
        # so go over each slot first and check that
        # it is empty
        total = self.num_filled + num_pmems
        for i in range(inputs.num_slots):
            if self.num_filled == total:
                break
            if self[i] is None:
                self[i] = Pmem(i, current_mev, inputs.num_geoms,
                               inputs.num_diheds)
                self[i].setup(major=False)

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
    def energies(self):
        """The set of energy values in the ring that are not None."""
        energies = []
        for pmem in self:
            if pmem:
                energies += [e for e in pmem.energies if e is not None]
        return energies

    @property
    def rmsds(self):
        """The set of RMSD values in the ring that are not None."""
        rmsds = []
        for pmem in self:
            if pmem:
                rmsds += [r[2] for r in pmem.rmsds if r[2] is not None]
        return rmsds

    @property
    def best_pmem(self):
        """Return the slot number and fitness value for the pmem with the highest fitness.

        Notes
        -----
        Does not update fitness values. Call update_all_fitness
        from the fitness module first.

        """
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
        """Get the mean fitness of the ring."""
        fit_vals = [self[i].fitness for i in self.occupied if self[i].fitness is not None]
        try:
            mean_fitness = mean(fit_vals)
        except StatisticsError:
            return None
        return mean_fitness

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
        if count == 0:
            return None
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
        if count == 0:
            return None
        return total / count

    @property
    def median_fitness(self):
        """Get the median fitness of the ring."""
        fit_vals = [self[i].fitness for i in self.occupied if self[i].fitness is not None]
        try:
            med_fitness = median(fit_vals)
        except StatisticsError:
            return None
        return med_fitness

    @property
    def median_energy(self):
        energy_vals = []
        for slot in self.occupied:
            for energy in self[slot].energies:
                # None means energy could not be calculated
                if energy is not None:
                    energy_vals.append(energy)
        try:
            med_energy = median(energy_vals)
        except StatisticsError:
            return None
        return med_energy

    @property
    def median_rmsd(self):
        rmsd_vals = []
        for slot in self.occupied:
            for rmsd in self[slot].rmsds:
                # None means one or two invalid geometries
                if rmsd[2] is not None:
                    rmsd_vals.append(rmsd[2])
        try:
            med_rmsd = median(rmsd_vals)
        except StatisticsError:
            return None
        return med_rmsd

    @property
    def stdev_fitness(self):
        """Get the standard deviation for fitness in the ring."""
        fit_vals = [self[i].fitness for i in self.occupied if self[i].fitness is not None]
        try:
            stdev_fitness = stdev(fit_vals)
        except StatisticsError:
            return None
        return stdev_fitness

    @property
    def stdev_energy(self):
        energy_vals = []
        for slot in self.occupied:
            for energy in self[slot].energies:
                # None means energy could not be calculated
                if energy is not None:
                    energy_vals.append(energy)
        try:
            stdev_energy = stdev(energy_vals)
        except StatisticsError:
            return None
        return stdev_energy

    @property
    def stdev_rmsd(self):
        rmsd_vals = []
        for slot in self.occupied:
            for rmsd in self[slot].rmsds:
                # None means one or two invalid geometries
                if rmsd[2] is not None:
                    rmsd_vals.append(rmsd[2])
        try:
            stdev_rmsds = stdev(rmsd_vals)
        except StatisticsError:
            return None
        return stdev_rmsds
