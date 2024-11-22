# One-dimensional examples

## Sod Shock Tube

The Sod Shock tube is a classical Riemann problem used to test the accuracy of computational methods, in particular the capacity of numerical methods to handle shocks.

The initial conditions are very simple for this problem: a contact discontinuity separating a gas with different pressure and density, and zero velocity everywhere. In the standard case the density and pressure on the left are unity, The density on the right side of the diaphragm is 0.125 and the pressure is 0.1. More details about this classic case can be obtained from Laney and Toro books.
The solution is plotted at t=0.2 The exact solution is calculated using [ToroExact](https://github.com/tahandy/ToroExact)


The example is located in the folder `numerics/riemmann`. After compiling using `make` to run type

```bash
$ ./main1d.gnu.MPI.ex inputs
```
It should run very fast (<1 secs, depending on machine)
The results can be seen by

```bash
$ python plot.py
```

(images/.png)
**Skew-symmetric** 4th order, Shock Capturing and damping with constants **C2=0.1** and **C4=0.016**  N=96
(images/num_rie_skew_damp.png)


**Skew-symmetric** 4th order,Shock Capturing and damping with constants **C2=0.5** and **C4=0.016**  N=96
(images/num_rie_skew_damp2.png)


**MUSCL** Riemann solver 2th order  N=96
(images/num_rie_muscl.png)


## Shu-Osher Problem


The problem corresponds to an  Mach 3 shock propagating into a field with small density disturbances, which can be interpreted as entropy disturbances. The solution is compared against a reference solution obtained using a fine mesh. 

The example is located in the folder `numerics/shu`. After compiling using `make` to run type

```bash
$ ./main1d.gnu.MPI.ex inputs_shu
```

Results can be seen using
```bash
$ python plot.py
```

The table belows shows an approximate relative cost (gcc version 11, 1 core, CFL fix 0.3),
in the same CPU

| Run                       | Mesh         | Run Time | 
| --------------------------| -------------|----------| 
| Reference MUSCL           | 8912         |  49.337  | 
| MUSCL                     | 256          |  0.098   | 
| Skew 4 *Cs=0.1 Cd=0.016*  | 256          |  0.08    | 

**MUSCL** Riemann solver 2th order  **N=256**
(images/num_shu_muscl.png)

**Skew**  solver 4th order  **N=256**, *Cs=0.1* *Cd=0.016*

(images/num_shu_skew.png)

**Skew**  solver 4th order  **N=512**, *Cs=0.1* *Cd=0.016*
(images/num_shu_skew2.png)