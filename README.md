SparseCovReg Analysis
================

This GitHub repository contains the codes for reproducing the analysis in ''High-dimensional covariance regression with application to co-expression QTL detection''.

## Code

This folder contains the codes for reproducing the analysis in Sections 4, 5 and the supplementary material of the paper (simulation studies and real data studies).

1.  simul.R: All results in the simulation study (Section 4) can be reproduced by this file.
2.  simul_supp.R: All results in the supplementary material (Section S9) can be reproduced by this file.
3.  simul_CV.R: All results in the simulation study (Section 4) can be reproduced by this file.

## How does the process work?

### Step 1

Author(s) can create a public GitHub repository in their own GitHub account
by using this template repository. This template contains a basic 
skeletal structure to help authors structure their code and analyses for their 
JASA publication. Creating a repository with the template can be done in the following way: 

Click on the "Use this template" button for [this GitHub template repository](https://github.com/jasa-acs/repro-template). (You'll need to be signed in to a GitHub account in order to see the button.)

![Click template button](https://docs.github.com/assets/cb-36544/images/help/repository/use-this-template-button.png)

From there, author(s) can [follow these instructions](https://docs.github.com/en/repositories/creating-and-managing-repositories/creating-a-repository-from-a-template). However do not optionally select "**Include all branches**" as you do not need this for your own projects. 


### Step 2

The author(s) can then directly edit (or replace) the manuscript template files in their own GitHub repository. Author(s) can also add their own data, code, and other files as needed. 

For guidance on getting started with git, we recommend the [Happy with git r](https://happygitwithr.com) tutorials.

**Importantly, the authors should provide an overview of how to carry
out the analyses presented in their manuscript in the `README.md` of their
repository, replacing the content in this file.** This overview would
generally refer to scripts/code files that execute the analyses and are
placed either in the main directory or the `/code` subdirectory. The
*Workflow* section of the ACC form should refer to this README.md as
containing the instructions for how to reproduce the analyses.

### Step 3

Author(s) use `git commit` to track changes over time and use `git push`
to push changes to a repository on the author(s) personal GitHub
account.

### Step 4

Author(s) submit a link to their GitHub repository as part of the [JASA
Reproducibility review process](https://jasa-acs.github.io/repro-guide/),
required upon submission of an invited revision.

### Step 5

JASA Associate Editors for Reproducibility will review the materials in
the GitHub repository of the authors and submit a
reproducibility review as part of the standard JASA review process.
Authors have the opportunity to respond to the review by making changes
and pushing their changes to their personal GitHub repository.

### Step 6

Once the manuscript is accepted, the materials in the author(s) personal
GitHub repository will be copied to the [JASA repository](https://github.com/jasa-acs).

## Reproducibility materials file structure

This template provides a suggested file structure for a JASA submission, but authors are free
to modify this structure.

The suggested components are as follows. Directories in the submission may have subdirectories to
further organize the materials.

1.  A `README.md` file - This file gives a short description of the
    paper and an overview of how to carry out the analyses presented in their manuscript.
2.  A `manuscript` directory - This directory will generally hold the source files
    (often LaTeX or Rmd) for the manuscript and any files directly related to the
    generation of the manuscript, including figure files.
3.  A `data` directory - This directory will generally hold the real data files 
    (or facsimile versions of them in place of confidential data) and simulated data files.
    See `data/README.md` for more details. 
4.  A `code` directory - This directory will generally hold 
    source code files that contain the core code to implement the method and various utility/auxiliary functions.
5.  An `output` directory - This directory will generally hold objects derived
    from computations, including results of simulations or real data analyses. See `output/README.md` for more details.

## Guidance on the use of reproducible environments

Submissions may include the use of reproducible environments capturing
state of a machine generating manuscript artifacts and even the
manuscript itself. Here we discuss two types of reproducible
environments and their use. Both virtual and package environments may be
put in the `code` directory.

### Package environments

Package environments capture the set of packages used by a programming
language needed to generate output. The R programming language has
`renv`, `switchr` and others to accomplish this, Python has `venv`,
`conda` and others, and Julia has native support (through the `Pkg`
package). When submitting these types of environments, the following are
suggested.

1.  Clearly indicate (in the overall `README.md`) the language(s) used (including version) 
    and the package environment tool used (e.g., `renv`, `conda`).
2.  Use a single package environment for all reproducible content.
3.  Prefer packages from package archives (CRAN, Bioconductor,
    RForge.net for example).
4.  If you use packages from a code repository (GitHub, GitLab, etc.)
    then use a release version if possible, or indicate the commit used. You could also consider
    forking the repository and providing a release.

### Virtual environments

Virtual environments such as Docker and Singlarity capture
the entire computing environment in which computations were performed.
In general, they are a more robust solution, capable of taking a
“snapshot” of a machine, including any system-level utilities and
external libraries needed to perform your computations. They have the
advantage that reproducing materials means running the virtual
environment, rather than recreating the programming language environment.
If using a virtual environment, we ask that 
you provide a definition file (e.g., a Dockerfile) or (perhaps better)
a link to an image in a standard online registry, such as DockerHub.

## References

Gentleman, Robert, and Duncan Temple Lang. “[Statistical Analyses and
Reproducible
Research](http://biostats.bepress.com/cgi/viewcontent.cgi?article=1001&context=bioconductor).”
(2004).

Gentleman, Robert. “[Reproducible research: a bioinformatics case
study](https://www.degruyter.com/document/doi/10.2202/1544-6115.1034/html).”
Statistical applications in genetics and molecular biology 4.1 (2005).

Marwick, Ben, and Bryan, Jennifer, and Attali, Dean, and Hollister,
Jeffrey W. [rrrpkg Github Page](https://github.com/ropensci/rrrpkg).
