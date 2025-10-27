![Example WMLES of a ramjet combustor. From top to bottom: temperature, pressure, density gradient, and streamwise velocity.](cerisse_ramjet.png)

# Cerisse
**A research computational solver for Large Eddy Simulation and Direct Numerical Simulation of compressible turbulent reactive flows**

Version 1 created by: Enson Un, Salvador Navarro-Martinez


## Capabilities
- Fifth-order low-dissipation WENO[^1] and TENO[^2] schemes
- LES/PDF model solved using the Eulerian Stochastic Fields approach[^3]
- Wall-modelled Stochastic Fields and wall-modelled LES with advanced SGS models
- Adaptive mesh refinement powered by AMReX[^4]
- Realistic thermodynamic and chemical properties using PelePhysics[^5]


## Getting the code
1. Clone the repository
    ```bash
    git clone git@github.com:salvadornm/cerisse.git
    cd cerisse
    ```
2. Checkout to the `cerisse1` branch
    ```bash    
    git checkout cerisse1
    ```
3. Get all submodules
    ```bash    
    git submodule init
    git submodule update
    ```
4. (Optional) Disable some PelePhysics constraints and aborts
    ```bash
    cd EB_CNS/Tools/scripts
    ./pp_modify.sh
    ```

Or download a release from [here](https://github.com/salvadornm/cerisse/releases).


## Compiling and running
We use GNU Make system for generating executables. The CMake system is not being maintained.
1. Go to an exec folder, for example,
    ```bash
    cd EB_CNS/Exec/ShockReflect
    ```
2. Compile the third-party SUNDIALS library first, which requires CMake. This only need to be done once per compiler settings combination.
    ```bash
    make TPL
    ```
3. Compile the executable. You may use the `-j` flag to build multiple jobs in parallel.
    ```bash
    make -j8
    ```
4. Run the executable, for example,
    ```bash
    mpirun -np 8 ./Cerisse2d.gnu.MPI.ex inputs
    ```

- See `prob_parm.H`, `prob.H`, and `prob.cpp` for the problem definition
- See `inputs` for runtime options. See [`default_parm.H`](EB_CNS/Source/default_parm.H) for more options and their default values.
- See `GNUmakefile` for compile-time options
- More information about the solver can be found in [`docs`](docs). To generate the docs locally, do `mkdocs serve` in the parent directory and point the browser to [127.0.0.1.8000](http://127.0.0.1:8000). You may also find [cerrisse2 docs](hslesdnsrf.gitbook.io/cerisse-docs) useful, but bear in mind that the two solvers are different.


## Citation
This open-source software is distributed under the [BSD3 license](LICENSE). If you use it in your work, please cite our [Combustion and Flame article](https://doi.org/10.1016/j.combustflame.2024.113897)

```
@article{Un2025,   
   author = {Tin-Hang Un and Salvador Navarro-Martinez},
   doi = {10.1016/j.combustflame.2024.113897},
   journal = {Combustion and Flame},
   month = {2},
   pages = {113897},
   title = {Stochastic fields with adaptive mesh refinement for high-speed turbulent combustion},
   volume = {272},
   year = {2025}
}
```

We recommend to cite as well AMReX and PelePhysics as in the references.


## References
[^1]: Castro *et al.* (2011). [High order weighted essentially non-oscillatory WENO-Z schemes for hyperbolic conservation laws](https://doi.org/10.1016/j.jcp.2010.11.028). *JCP*, 230(5), 1766–1792.

    Balsara *et al.* (2024). [Efficient Alternative Finite Difference WENO Schemes for Hyperbolic Conservation Laws](https://doi.org/10.1007/s42967-023-00360-z). *Comm. Appl. Math. Comput.*, 1–54. 
[^2]: Fu *et al.* (2016) [A family of high-order targeted ENO schemes for compressible-fluid simulations](https://doi.org/10.1016/j.jcp.2015.10.037). *JCP*, 305, 333–359.
[^3]: Un & Navarro-Martinez (2025). [Stochastic fields with adaptive mesh refinement for high-speed turbulent combustion](https://doi.org/10.1016/j.combustflame.2024.113897). *Combust. Flame*, 272, 113897.
[^4]: Zhang *et al.* (2019). [AMReX: a framework for block-structured adaptive mesh refinement](https://doi.org/10.21105/JOSS.01370). *JOSS*, 4(37), 1370.
[^5]: Henry de Frahan *et al.* (2024). [The Pele Simulation Suite for Reacting Flows at Exascale](https://doi.org/10.1137/1.9781611977967.2). *Proceedings of the 2024 SIAM Conference on Parallel Processing for Scientific Computing*, 13–25.
