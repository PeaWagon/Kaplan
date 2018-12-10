# Kaplan

Project for CAS741.

Developer Name: Jen Garner

This project implements conformer searching using evolutionary computation.

The main folders in the project are:
* docs (documentation)
* kaplan (source code)

There is also a test folder within the kaplan directory that contains:
* testfiles
* jupyter-notebooks

## Dependencies

Kaplan has the following dependencies. Note sublists indicate that Kaplan
does not directly import the dependency, but it is needed by its parent
in the list.

* [python 3.6](https://www.python.org/downloads/)
* [numpy](https://pypi.org/project/numpy/)
* [rmsd](https://github.com/charnley/rmsd)
    * [scipy](https://pypi.org/project/scipy/)
* [psi4]()
* [openbabel]()
* vetee
    * [pubchempy](https://pubchempy.readthedocs.io/en/latest/guide/install.html)

## How to install

### Installing Dependencies

The recommended installation process involves installing [conda or miniconda](https://conda.io/docs/user-guide/install/download.html#anaconda-or-miniconda).

Generate a new conda environment:  
`$ conda create -n kenv python=3.6 numpy`

Turn on the environment:  
`$ source activate kenv`

Your prompt should now reflect the environment name:  
`(kenv) $`

With the environment active, install the dependencies:  
1. psi4  
`(kenv) $ conda install -c psi4 psi4`  
2. openbabel  
`(kenv) $ conda install -c openbabel openbabel`  
3. pubchempy  
`(kenv) $ conda install -c mcs07 pubchempy`  
4. rmsd  
`(kenv) $ pip install rmsd`  
5. vetee  
`(kenv) $ pip install -i https://pypi.anaconda.org/kumrud/simple vetee`   

### Installing Kaplan

Do a git clone of the kaplan repo:  
`(kenv) $ git clone https://github.com/PeaWagon/Kaplan.git`

Go into the Kaplan directory:  
`(kenv) $ cd Kaplan`

Then do a pip install (development):  
`(kenv) $ pip install -e ./`

### Working with the environment 

To see installed packages:  
`(kenv) $ conda list`  

Turn off the environment:  
`(kenv) $ source deactivate kenv`  
`$`

### Uninstall Kaplan

`(kenv) $ pip uninstall kaplan`






