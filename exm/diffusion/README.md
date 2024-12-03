## Diffusion

This test examines AMR behavior in a temperature diffusion problem with an initial condition defined as :

$$
T =  T_0 + \Delta T  \exp(-r/\delta)^2
$$

with T0=300 and DT=100. $$delta$$ somehow determines the width
of initial pulse and periodic boundary conditions.

The set-up is based on HAMISH validation case

https://www.ukctrf.com/index.php/benchmarking-of-the-new-software/


Tested       |          grid | comment
:----------- |:-------------:| -----------:
gcc 14.2     | **64x64**(2 levels)        |   200 steps (t=0.02), AMR

This test case includes only a diffusion term and, as a result, does not conserve mass, momentum, or energy. The initial temperature pulse generates a non-periodic velocity field. This example is designed to assess the performance of the diffusion and heat flux modeling, rather than to represent a real-world scenario. It is important to note that the simulation will fail if heat is allowed to reach the boundaries.