# Cerisse Documentation

AMREX-based Navier-Stokes Reactive compressible LES code 

Created by:

Enson Un

Salvador Navarro-Martinez

version 0.0 November 2023


## Folder layout

    EB_CNS/         # main code
    Submodules/     # required modules: AMREX, PelePhysics
    docs/           # documentation folder


## Pre-requisites

1. **C++ compiler** 
A compiler with [C++17](https://en.wikipedia.org/wiki/C%2B%2B17) standard is required.
Examples include **gcc > 8**, **clang >7**  and more

2. **GNU Makefile** 
It will usually be installed by default in most systems (Mac OS/Linux)

3. **AMREX** 
AMR libraries [AMREX](https://amrex-codes.github.io/amrex/)
It will be installed if not present. Check Quickrun.

4. **PelePhysics**
Based on
[PeleC](https://github.com/AMReX-Combustion/PeleC)
It will be installed if not present. Check Quickrun.

5. **MPI libraries** 
(optional) Rrequired for parallel simulations
Similarly CUDA/OPenMPI may be required for more advanced parallelization strategies.

6. **cmake** 
(optional) reuired for some installation options (mostly related to GPU and chemistry). 
Easy to install,  version reuired **>3.2**

7. **Visualization**
Cerisse/AMREx/PeleC format is supported by
[VisIt](https://visit-dav.github.io/visit-website/), 
[Paraview](https://www.paraview.org), 
[yt](https://yt-project.org) (allows Python)
and check for more options
[AMReX Visualization](https://amrex-codes.github.io/amrex/docs_html/Visualization.html)

## Installation

Download the latest version on the master branch using GitHub CLI

```
$ gh repo clone salvadornm/cerisse
```

or Download the latest release from [GitHub](https://github.com/salvadornm/cerisse/releases)



### Documentation Editing
For help editing the documentation visit [mkdocs.org](https://www.mkdocs.org). To generate the docs locally `mkdocs serve`
and point the browser to [127.0.0.1.8000](http://127.0.0.1:8000)
You will need to install Math extension for Python-markdown

```
$ pip install python-markdown math
```

