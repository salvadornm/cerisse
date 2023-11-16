# Cerisse Documentation

AMREX-based Navier-Stokes Reactive LES code

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
Examples include gcc > 8, clang >7  and more

2. **MPI libraries** (only required for parallel simulations)

3. **GNU Makefile** 
It will usually be installed by default in most systems (Mac OS/Linux)

4. **cmake** 
Easy to install,  version >3,2


5. **AMREX** 
AMR libraries [AMREX](https://amrex-codes.github.io/amrex/)

6. **PelePhysics**
Based on
[PeleC](https://github.com/AMReX-Combustion/PeleC)

## Visualization 

The plotfile format supported are
[VisIt](https://visit-dav.github.io/visit-website/), 
[Paraview](https://www.paraview.org), 
[yt](https://yt-project.org) (allows Python)
and check for more options
[AMReX Visualization](https://amrex-codes.github.io/amrex/docs_html/Visualization.html)

## Quick Installation

Dowload the latest version on the master branch

Using GitHub CLI

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

