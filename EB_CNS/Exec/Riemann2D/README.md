# 2D Riemann Problem

This is a standard test case for Euler equation solvers. It can be used to test the reconstruction scheme, AMR, and numerical stability. The initial conditions are adopted from [Kurganov & Tadmor (2002)](https://doi.org/10.1002/num.10025).

The following figure shows the all 19 configurations in the reference. All cases are run by TENO5 scheme with characteristic variable interpolation to ensure minimal numerical dissipation.

![](all_configs.png)