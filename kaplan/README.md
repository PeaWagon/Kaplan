# Running Kaplan

Kaplan requires the user to select inputs. Most inputs have
a default setting and so do not need to be set. There are
two types of inputs:  molecular inputs and evolutionary algorithm
(EA) inputs.

## Basic Program Execution

1. Activate the environment where Kaplan is installed.
2. Import the run_kaplan function from the control module.
3. Write an inputs dictionary.
4. Run the run_kaplan function by passing it the inputs dictionary.

At minimum, the user must specify the name of the molecule for
which to perform a conformer search.

### Detailed Examples

Running the default program execution with propane:

```bash
conda activate kaplan_env
(kaplan_env) python
>>> from kaplan.control import run_kaplan
>>> input_dict = {"struct_input": "propane"}
>>> run_kaplan(input_dict)
```

To run Kaplan with a different input type (example: an xyz file), the
struct_type argument can be set:
`>>> input_dict = {"struct_input": "/home/user/1,3-butadiene.xyz", "struct_type": "xyz"}`
The other options for struct_type are specified by:

```python
>>> from kaplan.inputs import Inputs
>>> inputs = Inputs()
>>> print(inputs._avail_structs)
{'inchi', 'smiles', 'cid', 'xyz', 'inchikey', 'com', 'name'}
```

By default, Kaplan saves the final python ring object using pickle.
If you don't want to save any objects, set the save argument
as False like so:
`>>> run_kaplan(input_dict, save=False)`

The default location to save the output is in the current working
directory, under kaplan_output/job_x_molname, where x is a number
and molname is the name of the molecule (or input filename, CID,
InChIkey etc.). Each time kaplan is run from the same directory,
a new job directory is generated. To change the output location,
add a new argument to the input_dict (output_dir), as follows:

```python
>>> input_dict = {"struct_input": "propane", "output_dir": "/home/user/kaplan_jobs"}
>>> run_kaplan(input_dict)
```

The run_kaplan function also accepts a Kaplan inputs object as input
rather than the input_dict.

```python
>>> from kaplan.inputs import Inputs
>>> inputs = Inputs()
>>> inputs.update_inputs({"struct_input": "propane"})
>>> run_kaplan(inputs)
```

The user can also specify an inputs.pickle file. By default, the new output
is written to a new directory.

```python
>>> old_inputs = "/home/user/kaplan_jobs/kaplan_output/job_0_propane/inputs.pickle"
>>> run_kaplan(old_inputs)
```

To keep the same output directory, specify the new_dir argument as False:

```python
>>> old_inputs = "/home/user/kaplan_jobs/kaplan_output/job_0_propane/inputs.pickle"
>>> run_kaplan(old_inputs, new_dir=False)
```

The user can also specify an old ring object (if they wish to continue a job
with more mating events, for example):

```python
>>> old_inputs = "/home/user/kaplan_jobs/kaplan_output/job_0_propane/inputs.pickle"
>>> old_ring = "/home/user/kaplan_jobs/kaplan_output/job_0_propane/ring.pickle"
>>> run_kaplan(old_inputs, ring=old_ring)
```

## 1: EA Inputs

These are the inputs relating to the evolutionary algorithm. Defaults
are indicated in parentheses ().

* **num_mevs (100)**: number of mating events to do
* **num_slots (50)**: number of ring slots total
* **init_popsize (10)**: number of starting filled slots in the ring
* **num_geoms (5)**: number of conformers to find
* **num_muts (num_dihed * num_geoms // 5)**: the maximum number of mutations to apply to one solution instance
* **num_swaps (num_geoms // 2)**: the maximum number of swaps to perform between
two sets of geometries
* **num_cross (num_geoms // 2)**: the maximum number of crossovers
to perform between two sets of geometries
* **mating_rad (3)**: the maximum distance (in number of slots) that
a new population member can be placed away from its parent. Also decides
the size of the tournament.
* **fit_form (0)**: only fit_form 0 is supported at the moment
* **coef_energy (0.5)**: the coefficient of the energy summation in
the fitness function
* **coef_rmsd (0.5)**: the coefficient of the root-mean-square
deviation summation in the fitness function
* **asteroid (0.0)**: the percent chance per mating event that the
asteroid operator will be applied to the ring
* **plague (0.0)**: the percent chance per mating event that the
plague operator will be applied to the ring
* **agathic (0.0)**: the percent chance per mating event that the
agathic operator will be applied to the ring
* **deluge (0.0)**: the percent chance per mating event that the
deluge operator will be applied to the ring

## 2: Molecular Inputs

These are the inputs relating to the genetic algorithm. Defaults
are indicated in parentheses (). Required arguments are specified
using a star *.

* **qcm (hf for psi4, mmff94 for openbabel)**: quantum chemical method or forcefield used for energy calculations
* **basis (sto-3g)**: basis set to use
* **struct_input\***: should be a file name, string, or cid
* **struct_type\***: one of xyz, com, glog, smiles, cid, or name
* **charge\***: charge of the input molecule
* **multip\***: multiplicity of the input molecule

Other molecular attributes:

* **prog (psi4)**: the program to run the energy calculation. Can
be psi4 or openbabel

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
