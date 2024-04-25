# Cerisse++



--8<-- "README.md"


## Pre-requisites

1. **C++ compiler** 
A compiler with [C++17](https://en.wikipedia.org/wiki/C%2B%2B17) standard is required.
Examples include **gcc > 8**, **clang >7**  and more

2. **GNU Makefile** 
It will usually be installed by default in most systems (Mac OS/Linux)

3. **MPI libraries** 
(optional) required for parallel simulations.
Similarly CUDA/OPenMPI may be required for more advanced parallelization strategies.

4. **cmake** 
(optional) required for some installation options (mostly related to GPU and chemistry). 
Easy to install, version required **>3.2**

5. **AMReX** 
AMR libraries [AMREX](https://amrex-codes.github.io/amrex/)
It will be installed if not present. Check [Quickrun](quickrun.md) for installation within Cerisse (recommended)

6. **PelePhysics**
Based on [PeleC](https://github.com/AMReX-Combustion/PeleC)
It will be installed if not present. Check [Quickrun](quickrun.md) for installation within Cerisse (recommended)

7. **Visualization**
Cerisse/AMREx/PeleC format is supported by
[VisIt](https://visit-dav.github.io/visit-website/), 
[Paraview](https://www.paraview.org), 
[yt](https://yt-project.org) (allows Python)
and check for more options
[AMReX Visualization](https://amrex-codes.github.io/amrex/docs_html/Visualization.html)


## Start

To finih installation see [QuickStart](quickrun.md)

## Folder layout

    cerisse/           # main code
    ├───exm/           # simple example cases
    ├───tst/           # test cases
    ├───src/           # souce code
    └───bin/           # useful scripts
    lib/               # required modules
    ├───amrex/         # AMR data structure, backends
    ├───PelePhysics    # EoS, transport, reaction
    └───sundials/      # chemical integrators
    docs/           # documentation folder




