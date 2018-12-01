# Kaplan
Project for CAS741.

Developer Name: Jen Garner

This project implements conformer searching using evolutionary computation.

There are four folders in the project:
* docs (documentation)
* refs (references/literature)
* src (source code)
* test (test files)


## How to install

Generate a new conda environment:
$ conda create -n kenv python=3.6 numpy

Turn on the environment:
$ source activate kenv
(kenv) $

Turn off the environment:
(kenv) $ source deactivate kenv
$

### Other Dependencies

With the environment active:
1. psi4
(kenv) $ conda install -c psi4 psi4
2. openbabel
(kenv) $ conda install -c openbabel openbabel
3. pubchempy
(kenv) $ conda install -c mcs07 pubchempy

