
These tests (located in ```exm/numerics```) consist of two classical tests:
the convection of a vortex (same configuration a sin ```exm/covo```)
and a 1D Riemann problem. The test are good to test behaviour 
of numerical schemes, regarding shock capturing and dissipation of a vortex.


The problem data is defined in two folders


### Convective-Vortex

**Skew-symmetric** 4th order no shock capturing nor damping  N=128 x 128
![covnum](images/num_covo_skew_nodamp.png)


**Skew-symmetric** 4th order,Shock Capturing and high freq damping with constants **C=0.5** and **C=0.0016**  N=128 x 128
![covnum2](images/num_covo_skew_damp.png)


**MUSCL** Riemann solver 2th order  N=128 x 128
![covonum3](images/num_covo_muscl.png)


### Riemann Problem


**Skew-symmetric** 4th order,Shock Capturing and high freq damping with constants **C=0.5** and **C=0.0016**  N=96
![rienum1](images/num_rie_skew_damp.png)


**Skew-symmetric** 4th order,Shock Capturing and high freq damping with constants **C=0.1** and **C=0.00016**  N=96
![rienum2](images/num_rie_skew_damp2.png)


**MUSCL** Riemann solver 2th order  N=96
![rienum3](images/num_rie_muscl.png)


