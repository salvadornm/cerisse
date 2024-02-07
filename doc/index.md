# Cerisse 


--8<-- "README.md"

## Folder layout

    EB_CNS/         # main code
    ├───Examples/      # simple example cases
    ├───Exec/          # more advanced test cases
    ├───Source/        # source code
    └───Tools/         # useful tools
    Submodules/     # required modules
    ├───amrex/         # AMR data structure, backends
    ├───PelePhysics    # EoS, transport, reaction
    └───sundials/      # chemical integrators
    docs/           # documentation folder

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

## Installation

Download the latest version on the master branch using GitHub CLI

```bash
$ gh repo clone salvadornm/cerisse
```

or Download the latest release from [GitHub](https://github.com/salvadornm/cerisse/releases)

## Documentation Editing

For help editing the documentation visit [mkdocs.org](https://www.mkdocs.org). To generate the docs locally type in the parent directory: `mkdocs serve`
and point the browser to [127.0.0.1.8000](http://127.0.0.1:8000)

You will need to install the `python-markdown-math` extension for rendering equations and the `markdown-callouts` extension for correctly displaying the warning and note blocks. All requirements can be installed automatically using

```bash
$ pip install -r docs/requirements.txt
```

You may need to install

```bash
$ pip install pip-tools
```

if you add new markdown extensions, edit the `requirements.in`  file under `docs/`

```bash
$ pip-compile requirements.in
```

