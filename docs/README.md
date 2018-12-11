# Documentation

This directory contains subdirectories for Kaplan's 
documentation. The templates for all of the documents were
taken from a gitlab repository -
[cas741](https://gitlab.cas.mcmaster.ca/smiths/cas741), authored
by [Dr. Spencer Smith](https://github.com/smiths). Another example of
this style of documentation can be found
[here](https://github.com/smiths/swhs).

### How to make the documentation

There are pdfs available for the documentation in the subdirectories.
If changes are made to the tex files, then you will need to regenerate the pdfs.
In order to generate the pdfs, you need to have latex (texlive), make, bibtex,
and pdflatex installed. Then, run the following command:

`(kenv) $ make -f Makefile`

An IDE for editing tex files is recommended (such as texstudio or texmaker).

To get rid of all of the extra files, run the clean command:

`(kenv) $ make -f Makefile clean`

### Design

This directory contains the Module Instance Specification (MIS)
document and the Module Guide (MG). Together, these documents
outline how the code is pieced together and provide specific
inputs and outputs for Kaplan modules.

### ProblemStatement

This directory contains a general abstract overview of the
problem that Kaplan solves.

### SRS

This directory contains the Software Requirements Specification
(SRS) document, which explains the purpose of Kaplan in abstract
terms. The SRS provides the theory, symbols, units, formulas,
and data formats necessary to answer the problem statement. The
requirements (functional and non-functional) given in this document
are intended to be addressed later on in the documentation (such as
in the design and in the testing).

### VnVPlan

This directory is for the Verification & Validation (VnV) Plan,
and it contains two main pieces of documentation:
1. SysVnVPlan addresses the system as a whole and follows the
SRS. Here we don't assume that it is known exactly how the
system will work, we just try to test its outputs and ensure
that the correct inputs are being given such that the
requirements can be satisfied.
2. UnitVnVPlan addresses the modules and follows the MIS. Here
we know exactly how the program works, and we want to write
tests to catch problems that might arise from implementation.

### VnVReport

This directory contains two reports - one for the SystVnVPlan
and the other for the UnitVnVPlan. It is used to document
progress that has been made in ensuring Kaplan's requirements
are being satisfied.

