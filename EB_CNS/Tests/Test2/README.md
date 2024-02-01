The Shu-Osher Problem is a one-dimensional shock-turbulence interaction in which a shock propagates into a density field with artificial fluctuations. The goal is to test the capability to accurately capture a shock wave, its interaction with an unsteady density field, and the waves propagating downstream of the shock.

The initial conditions on a domain [0,10] are:


$$
(\rho,u,p)  = \left\{ 
    \begin{array}{cl}
      (3.857143, 2.629369, 10.3333 ) & \quad \textrm{ if } x \geq 1 \\
      (1 + 0.2 \sin (5 x), 1 ) & \quad \textrm{ otherwise}
    \end{array}  
    \right.
$$

After compiling using `make` to run type
```
$ ./Cerisse1d.gnu.ex inputs
```

The figure shows three resolutiona

![test1plot](images/Test4_compare_grid.png)

**NOT FULLY TESTED**


