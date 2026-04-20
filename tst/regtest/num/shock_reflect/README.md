
This test ins a Mach 10 shock wave incident at 30 degrees (with adiabatic rati of 1.4).
 The initial conditions are given in [Berger & Colella (1989)](https://doi.org/10.1016/0021-9991(89)90035-1). The top boundary condition varies with time to follow the shock following the
unshocked speed. This test uses a tilted domain without EB.

Tested       |          grid | comment
:----------- |:-------------:| -----------:
gcc 11.4     | **56x16**(3 levels)        |  Euler,MUSCL 200 teps (t=0.2), AMR



