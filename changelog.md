# New in Version 1.3a.1

* remove unnecessary import statements from pmem, extinction, and tools modules.
* update fitness function for absolute fitness such that invalid geometries are penalised if the energies are positive for that pmem.
* fix issue with write_coords function where 7e-06 appears instead of 0.000007 in xyz files (thus breaking rmsd package)
* add function in optimise module to optimise Openbabel molecule object with RDKit (implementation to follow in later version)
* add function to geometry module to produce a set of minimum dihedral angles (one per rotatable bond), using the full set of dihedral angles from the Openbabel molecule object as a starting point
* remove GOpt dependency

## New in Version 1.3a.0

* Default `prog` changed to Openbabel from Psi4
* added new input option, `stop_at_conv`, which checks the best `pmem` for improvement every `stop_at_conv` mating events
  * this option ignores the `num_mevs` input option
  * if there is no change to the energies or RMSDs of the best pmem after `stop_at_conv` mating events, then the program terminates
* add mating event number to the last stats file update
* if flat molecule is detected, coordinates are generated from a SMILES string instead of InChI string, since some larger molecules were giving garbage coordinates when read as an InChI string by Openbabel.
  * Openbabel coordinates are scrambled for nonadecane (19 carbon chain)
  * Pubchem does not generate coordinates for this molecule (hence why it ends up being flat)
  * reading the InChI string with pybel results in completely garbage coordinates, but the SMILES string is okay
  * here is InChI from pubchem website for nonadecane: InChI=1S/C19H40/c1-3-5-7-9-11-13-15-17-19-18-16-14-12-10-8-6-4-2/h3-19H2,1-2H3
* make `name` a full `inputs` attribute (not just a property) to reduce times that Pubchem is queried during output plot generation
* fix bug where crossover fails because only one dihedral is present in the input molecule (cut index cannot be generated)
  * default behaviour is to check if `num_diheds` is 1, and (if true) set `num_cross` to zero
* add `n_point_crossover` function to the mutations module
* add input parameter `max_cross_points`, which represents how many cuts to make in `n_point_crossover`
  * the default for this parameter is the `num_diheds//3` (unless there are fewer than 3 dihedrals - for 2 dihedrals `max_cross_points` = 1 and for 1 dihedal crossover is turned off completely)
  * `n_point_crossover` is now used instead of the single point crossover that was previously implemented
* add input parameter `cross_points`, which is a list of indices, each in `range(1, num_diheds-1)`
  * the user can now select where exactly to perform crossover, otherwise the locations are chosen randomly from the above range
  * this change is meant to be used to keep important dihedrals (those that depend on one another) together during crossover
* add geometry test to make sure that GetX(), GetY(), and GetZ() return the same values as x(), y(), and z(), and to make sure sdf files (when turned into an obmol object) have the same geometry as get_coords
* add geometry optimisation in optimise module. Available for Openbabel and Psi4. Also added tests for this function to compare its results with the Openbabel binary, `obminimize`.
  * there are two types of optimise: major and minor
  * major optimisation occurs during the initial molecule setup (which can be turned off using the `opt_init_geom` input parameter - `True` by default)
  * major optimisation is also performed on the final best pmem in the ring during the call to the output module
  * minor optimisation is performed for each new pmem using a forcefield - the energy evaluation is performed based on the chosen program (Psi4 or Openbabel). This optimisation should get rid of overlapping atoms and issues where structures are invalid and cannot be read from Openbabel as xyz files.
  * minor optimisation is controlled via the input parameters `minor_tolerance`, `minor_maxsteps`, and `minor_sampling`. These inputs control the largest energy difference before convergence is achieved, the maximum steps to use of conjugate gradients before the optimisation is concluded, and how often (in terms of number of steps) to measure the change in energy, respectively. These same parameters are given for major optimisation, but these values only change the optimisation performed for Openbabel (since major optimisation for Psi4 is done using their optking module).
* add input parameter `avail_diheds`, which is an array of values from the half [-pi, pi). By default, numpy's linspace function is used to generate 16 values in this range. The user may add their own array as input such that only certain dihedral angles are chosen.
* fitness-related functions have been put in the fitness module and removed from the pmem and ring modules

## Possible Future Work

* add an input where the user can choose how to generate the initial coordinates (via Openbabel or by querying Pubchem, etc.) - this way, the error handling and input generation will be easier to follow
  * also, it may be a nice way to avoid the user needing an internet connection if they choose not to use Pubchem
  * cleanup of the inputs module is desired
* cut-down on dependencies on private repositories
* enable RMSD reorder method for final pmem
* ensure statistical relevance of data in the dihedrals heatmap; as a first start, shift dihedral angles such that the standard deviation in dihedral angle is minimised
* if there is a new release of Openbabel, re-write optimise function such that the changes to forcefield.cpp are accounted for (i.e. conjugate gradients should converge on its own - no need to calculate the energy delta)
* bug fixes are expected - a lot of changes were made to the code
* apply kabsch rotation matrix to final best pmem conformers so that their geometries are aligned with conf0.xyz
* keep track of whether an attempt has been made to calculate pmem energies/rmsds by initialising these attributes to None instead of empty lists; also don't need to store conformer indices in the RMSD (can make this a property), since pairs are reproducible using all_pairs_gen
* add the ability to make ramachandran plots to the tools module.
* counter in inputs to keep track of how many optimisations converged prior to maxsteps versus how many optimisations were run for maxsteps. This tracker would make it more obvious to the user if the calculations needed more time or if maxsteps was sufficient (tracker should be for major and minor).
* Paul suggested that the np.linspace function be changed to 12, 24, 36 steps such that 2pi/3 and pi/2 are included in the options.
