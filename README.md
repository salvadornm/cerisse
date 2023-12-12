# CERISSE
**A high-order adaptive mesh refinement solver for compressible turbulent reactive flows**

## Getting started
1. Clone the repository with all submodules
    ```bash
    git clone --recursive git@github.com:salvadornm/cerisse.git
    ```
3. Build the executable using Make
    ```bash
    cd EB_CNS/Exec/Riemann2D
    make
    ```
    You can speed up the process using the `-j` flag to specify the number of jobs to be run in parallel, e.g.
    ```bash
    make -j4
    ```
3. Run the executable (the executable name depends on the build configuration and system environment, you can find more details in the `GNUmakefile`)
    ```bash
    ./Cerisse2d.xxx.xxx.ex inputs
    ```
    or in parallel using MPI,
    ```bash
    mpiexec -np 8 ./Cerisse2d.xxx.xxx.ex inputs
    ```

## TODO List
- [] NSCBC (See aerofoil example, especially at high angle of attacks where vortex is shedded downstream)
- [] BC on EB
- [] LES: turbulent bulk viscosuty