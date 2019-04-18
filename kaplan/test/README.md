
***Kaplan Test Directory***

This directory is for testing purposes, and it contains
testing files for each module, a testfiles subdirectory,
and a jupyter-notebooks subdirectory.

**Testing Kaplan**

This directory contains the tests for Kaplan.
Each module has its own unit testing, under
test_modulename.py.

The testing is done via nosetests:

```
$ source activate kenv
(kenv) $ conda install nose
(kenv) $ nosetests kaplan
```

To test a single module, do:

`(kenv) $ nosetests kaplan.test.test_modulename`

For example, to test the mutations module, do:

`(kenv) $ nosetests kaplan.test.test_mutations`

To test a single test within a module, do:

`(kenv) $ nosetests kaplan.test.test_modulename:test_testname`

For example, to test the run_kaplan function in the gac module, do:

`(kenv) $ nosetests kaplan.test.test_gac:test_run_kaplan`

To run the tests without hiding printing:

`(kenv) $ nosetests --nocapture kaplan`

**Running jupyter notebooks**

To use the jupyter notebooks, you need to have jupyter installed.
Activate the conda environment where the other prerequisites are installed,
install jupyter, and then open the notebook:

```
$ source activate kenv
(kenv) $ conda install jupyter
(kenv) $ jupyter notebook
```
