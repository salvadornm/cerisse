
This test consist of the convection of a two-dimensional vortex in an ideal gas under M = 0.1 using Euler equations. The physical parameters are defined following Onera’s original definition

The vortex in the centre of the domain is defined by the stream function

$$
\Psi =  \Gamma \exp(-r/R)^2
$$

where *r* is relative to the initial positio of the vortex (x0,y0).


The velocity is 

$$
u = U_0 + \frac{\partial \Psi}{\partial y}   \; \; \; \; \; \; \; \;  v = -\frac{\partial \Psi}{\partial x} 
$$

and the pressure
$$
p =  p_0 - \frac{2 \Gamma^2 }{ R^2} \exp(-r/R)^2
$$

The solution is the actual vortex propagated in *x* and the solution repeats itself every flow-through time  L/U0.

The analytical solution is then (by differentiating the stream function):

$$
u = U_0  -  \frac{2 \Gamma }{ R^2} \exp(-r/R)^2 (y-y_0)
$$

and

$$
v =   \frac{2 \Gamma }{ R^2} \exp(-r/R)^2 (x-x_0)
$$


The example is in `exm/covo`

tested       |      grid     | comment
:----------- |:-------------:| -----------:
gcc 13.x (Mac)       | **128 x 128**        |  Euler, CFL=0.5, no-AMR, MPI

The results can be seen by

```
$ python plot.py
```

The results after five flow through times are

![covoplot](images/covo.png)

