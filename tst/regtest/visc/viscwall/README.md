
This is test is a 2-D laminar periodic channel flow.
The set-up is inspired on HAMISH validation case

https://www.ukctrf.com/index.php/benchmarking-of-the-new-software/


Tested       |          grid | comment
:----------- |:-------------:| -----------:
gcc 14.2     | **8x64**(2 levels)        | 


For plane Poiseuille flow (laminar flow between two parallel plates), the velocity profile is parabolic, and the maximum velocity is 1.5 the bulk velocity.
​The pressure gradient is a function of flow rate and bulk Reynolds number

$$
\frac{dp}{dx} = \frac{12 \rho Q^2}{L^3 \mbox{Re}_{bulk}}
$$

For a Reynolds number of 50, the pressure gradient is 0.24. The bulk velocity, density, and channel height are all set to unity, effectively using non-dimensional units.
DIFF number is 0.4 and bulk CFL is 0.32

The flow is initialized with a uniform bulk velocity of 1 and allowed to evolve until a parabolic profile is obtained. To compare it with the analytical solution, run the command: ```python ./plot.py``` 

![viscwall](images/viscwall.png)






