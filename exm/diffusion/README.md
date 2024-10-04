
This test ins a Mach 10 shock wave incident at 30 degrees (with adiabatic rati of 1.4).
 The initial conditions are given in [Berger & Colella (1989)](https://doi.org/10.1016/0021-9991(89)90035-1). The top boundary condition varies with time to follow the shock following the
unshocked speed. This test uses a tilted domain without EB.

Tested       |          grid | comment
:----------- |:-------------:| -----------:
gcc 11.4     | **56x16**(3 levels)        |  Euler,MUSCL 200 teps (t=0.2), AMR


### Boundary Condition Implementation

This test case is a godd example to incorporate complex boundary conditions using the ```bcnormal``` function

For example, to define manually the south boundary condition
in `prob.h`

```cpp
    case  2:  // SOUTH
     if (x[0] < 1. / 6.) {
      // post-shock conditions
      s_ext[URHO] = prob_parm.rho_l;
      s_ext[UMX]  = prob_parm.rho_l * prob_parm.u_l;
      s_ext[UMY]  = prob_parm.rho_l * prob_parm.v_l;
      s_ext[UMZ]  = 0.0;
      s_ext[UET]  =   prob_parm.eint_l + 0.5 *prob_parm.rho_l*
      (prob_parm.u_l * prob_parm.u_l + prob_parm.v_l * prob_parm.v_l);
      }
      else {
      // slip wall
      s_ext[URHO] = s_int[URHO];
      s_ext[UMX]  = s_int[UMX];
      s_ext[UMY]  = -s_int[UMY];
      s_ext[UMZ]  = s_int[UMZ];
      s_ext[UET]  = s_int[UET];
      }
      break;
```cpp

### Results
The density iso-contours and mesh refinment looks like:

![shock](images/shock_reflec.png)

Comparing 3 qand 6 refinement levels. Detail
 Due to the increased resolution, a  Kelvin-Helmholtz rollup can be seen along the principal slip line

![shock3](images/shock3.png)
![shock6](images/shock6.png)



