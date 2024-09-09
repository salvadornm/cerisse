(in progress)

A Mach 10 shock wave incident at 30 degrees. The initial conditions are given in [Berger & Colella (1989)](https://doi.org/10.1016/0021-9991(89)90035-1). The top boundary condition varies with time to follow the shock. The shock speed is given by 
$$
u_s = M c_s
$$
where $$M$$ and $$c_s$$ are Mach number and unshocked speed of sound, respectively. This version uses a tilted domain without EB.

Tested       |          grid | comment
:----------- |:-------------:| -----------:
gcc 12.2 Mac       | **56x16**(3 levels)        |  Euler, ??,  222 steps (t=0.2), AMR,  2 cores


### Boundary Condition Implementation

This test case is a godd example to incorporate complex boundary conditions using the ```bcnormal``` function




### Results
In Visit, the density iso-contours looks like:
