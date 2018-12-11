# Running Kaplan

Kaplan requires the user to write two input files:  
1. Molecular (mol) input file.  
2. Genetic algorithm (ga) input file.  
More help on how to write these files is given below.  

Then, the user must call the kaplan program using the gac
(genetic algorithm control) module with python
and provide the two input files as arguments, as so:

`(kenv) $ python gac.py mol_input_file.txt ga_input_file.txt`

## 1: mol input file

The molecular input file has the following parameters that must
be given:  
* **qcm**: quantum chemical method used for energy calculations
* **basis**: basis set to use
* **struct_input**: should be a file name, string, or cid
* **struct_type**: one of xyz, com, glog, smiles, cid, or name
* **prog**: only psi4 is supported at the moment
* **charge**: charge of the input molecule
* **multip**: multiplicity of the input molecule

## 2: ga input file

The genetic algorithm input file has the following parameters
that must be given:  
* **num_mevs**: number of mating events to do
* **num_slots**: number of ring slots total
* **num_filled**: number of starting filled slots in the ring
* **num_geoms**: number of conformers to find
* **num_atoms**: number of atoms for the input molecule
* **t_size**: the tournament size (how many slots to choose
from the ring during a mating event)
* **num_muts**: the maximum number of mutations to do on one
list of dihedral angles
* **num_swaps**: the maximum number of swaps to perform between
two sets of geometries
* **pmem_dist**: the maximum distance (in number of slots) that
a new population member can be placed away from its parent
* **fit_form**: only fit_form 0 is supported at the moment
* **coef_energy**: the coefficient of the energy summation in
the fitness function
* **coef_rmsd**: the coefficient of the root-mean-square
deviation summation in the fitness function

## Finding the output

The output is written to the kaplan/kaplan_output directory
by default. If you wish to change this location, open the
output.py file in the kaplan directory, and change the
variable OUTPUT_DIR to the desired location. Make sure
your version of OUTPUT_DIR appears above the run_output
function (and is the only defintion of OUTPUT_DIR).

## If you get an error

Here are some basic checks to do:  
* All of the dependencies should be installed (see main
repo page for instructions). The python version *must* be
at least 3.6 (otherwise you will get syntax errors for the
f-strings).
* If using an environment, the environment should be active.
* Example input files are available in the
kaplan/test/testfiles directory. Check that these inputs are
followed.
* Make sure that the order is molecule followed by genetic
algorithm input, otherwise the program will not run
successfully.
* Make sure the files are saved in the same directory, and
that the program is run from that same directory (otherwise
Kaplan won't be able to find the files).
* Raise a github issue on the Kaplan repo if you cannot get
the files to work.

