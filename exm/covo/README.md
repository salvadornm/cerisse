
This test (located in ```exm/covo```) consist of the convection of a two-dimensional vortex in an ideal gas under Mach = 0.1 using Euler equations. The physical parameters are defined following Onera’s original definition.

The vortex in the centre of the domain is defined by the stream function

$$
\Psi =  \Gamma \exp(-r/R)^2
$$

where *r* is relative to the initial positio of the vortex (x0,y0).
The velocity is therefore

$$
u = U_0 + \frac{\partial \Psi}{\partial y}   \; \; \; \; \; \; \; \;  v = -\frac{\partial \Psi}{\partial x} 
$$

and the pressure
$$
p =  p_0 - \frac{2 \Gamma^2 }{ R^2} \exp(-r/R)^2
$$

The solution is the actual vortex propagated in *x* and the solution repeats itself every flow-through time  L/U0.

The analytical solution can be obtained by differentiating the stream function:

$$
u = U_0  -  \frac{2 \Gamma }{ R^2} \exp(-r/R)^2 (y-y_0)
$$

and

$$
v =   \frac{2 \Gamma }{ R^2} \exp(-r/R)^2 (x-x_0)
$$


The parameters for the test are:

|            | value      | 
|------------|------------|
| p0         | 101300  Pa |
| rho0       |  1.1717047 kg/m3    |
| U0         | 35  m/s    |
| L          | 0.3112  m  |
| R          | L/20       |
| beta       | 0.04  |
|------------|------------|

where **beta* represents the strength of the vortex such that:

$$
\Gamma =  \beta U_0 R \sqrt{ \exp(1)}
$$






tested       |      grid     | comment
:----------- |:-------------:| -----------:
gcc 13.x (Mac)       | **128 x 128**        |  Euler, CFL=0.5, no-AMR, MPI

To run (with 2 MPI ranks)
```
$ mpirun -np 2 ./main2d.gnu.MPI.ex inputs
```


The results and comparison with theory can be seen by

```
$ python plot.py
```

Two prob.h files exist in the folder, the second order MUSCL Riemann solver (```prob.h_riemann```) and the skew-symmetric central differences (```prob.h_riemann```). Copy to ```prob.h``` to run the relevant scheme

The results after **five** flow through times for the Riemann-solver with MUSCL 2nd order reconstruction are

![covoplot](images/covo.png)

which show strong disspation.
Using the skew-symmetric 4th-order central scheme (without any high-frequency damping), the results exhibit less dissipation.

![covoplot2](images/covo_skew.png)

The distribution of vertical velocity can be seeen
![covot0](images/covo_vel.png)



