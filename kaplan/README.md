# Running Kaplan

Kaplan requires the user to select inputs. Most inputs have
a default setting and so do not need to be set. There are
two types of inputs:  molecular inputs and genetic algorithm
(GA) inputs.

### Basic Program Execution

First, open the terminal and change directory into where the
output files should be saved. At minimum, the user must specify
the molecule (via SMILES, name, cid, xyz file, or com file), its
charge, and its multiplicity. The example given below is for
the molecule propane, specified by name.
```
(kenv) $ python
>>> from kaplan.control import run_kaplan
>>> input_dict = {
        "struct_input": "propane",
        "struct_type": "name",
        `"charge": 0,
        "multip": 1,
        }  
>>> run_kaplan(input_dict)
```
By default, Kaplan saves the final python ring object using pickle.
If you don't want to save any objects, input the second argument
as False like so:
 
`>>> run_kaplan(input_dict, False)`  

The output is written to the present/current working directory
under kaplan_output/job_0 (for the first job). To specify a
specific output directory where kaplan_output/job_0 is generated,
enter that as a third argument:

`>>> run_kaplan(input_dict, True, "/home/user/mydir")`  


## 1: GA Inputs

These are the inputs relating to the genetic algorithm. Defaults
are indicated in parentheses ().  
* **num_mevs (10,000)**: number of mating events to do
* **num_slots (100)**: number of ring slots total
* **init_popsize (10)**: number of starting filled slots in the ring
* **num_geoms (5)**: number of conformers to find
* **num_muts (num_atoms/3)**: the maximum number of mutations to do on one
list of dihedral angles
* **num_swaps (num_geoms)**: the maximum number of swaps to perform between
two sets of geometries
* **mating_rad (3)**: the maximum distance (in number of slots) that
a new population member can be placed away from its parent. Also decides
the size of the tournament.
* **fit_form (0)**: only fit_form 0 is supported at the moment
* **coef_energy (0.5)**: the coefficient of the energy summation in
the fitness function
* **coef_rmsd (0.5)**: the coefficient of the root-mean-square
deviation summation in the fitness function

## 2: Molecular Inputs

These are the inputs relating to the genetic algorithm. Defaults
are indicated in parentheses (). Required arguments are specified
using a star *. 
* **qcm (hf)**: quantum chemical method used for energy calculations
* **basis (sto-3g)**: basis set to use
* **struct_input\***: should be a file name, string, or cid
* **struct_type\***: one of xyz, com, glog, smiles, cid, or name
* **charge\***: charge of the input molecule
* **multip\***: multiplicity of the input molecule

Other molecular attributes:
* **num_atoms**: number of atoms for the input molecule
* **prog**: only psi4 is supported at the moment

## If you get an error

Here are some basic checks to do:  
* All of the dependencies should be installed (see main
repo page for instructions). The python version *must* be
at least 3.7 (otherwise you will get syntax errors for the
f-strings and import errors).
* If using an environment, the environment should be active.
* Example input files are available in the
kaplan/test/testfiles directory.
* Raise a github issue on the Kaplan repo if you cannot get
the program to work.

