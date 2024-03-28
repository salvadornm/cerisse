
The Sod Shock tube is a classical Riemann problem used to test the accuracy of computational methods

The initial conditions are very simple for this problem: a contact discontinuity separating gas with different pressure and density, and zero velocity everywhere. In the standard case the density and pressure on the left are unity, The density on the right side of the diaphragm is 0.125 and the pressure is 0.1. More details about this classic case can be obtained from Laney[^1]  and Toro books.
The solution is plotted at t=0.2

tested       |      grid     | comment
:----------- |:-------------:| -----------:
gcc 11.3(Linux), 13.x (Mac)       | **200**        |  Euler, 200 steps, no-AMR

The exact solution is calculated using [ToroExact](https://github.com/tahandy/ToroExact)


After compiling using `make` to run type
```
$ ./main1d.gnu.MPI.ex inputs
```
It should run very fast (<1 secs, depending on machine)
The results can be seen by

```
$ python plt.py
```

![test1plot](images/test1.png)



[^1]: C.B. Laney *Computational Gas Dynamics*, Cambridge University Press, 1998.
[^2]: E.F. Toro. *Riemann solvers and numerical methods for fluid dynamics: a practical introduction*. Springer, Berlin, New York, 2009.

