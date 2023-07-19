# Supersonic turbulent flow in a plane channel with isothermal walls

## Setups
| Parameter            | Value                |
| -------------------- | -------------------- |
| Mach number $M = \frac{\langle u \rangle}{c_w}$ | 1.5 |
| Reynolds number $Re = \frac{\langle \rho \rangle \langle u \rangle H}{\mu_w}$ | 3000 |
| Prandtl number $Pr$  | 0.7 |
| Specific heat capacity ratio $\gamma$ | 1.4 |
| Dynamic viscosity $\mu = T^{0.7}$ | / |

Subscipts $w$ denote quantities at walls. $H=1$ is the channel half-height. The computational domain is $4\pi H \times 2H \times \frac{4\pi}{3}H$.

The flow is initialised with a laminar parabolic velocity profile $u = 1.5(H - y^2)$, $v = w = 0$, $T = \rho = 1$, with random velocity super-imposed.


## Analysis
Coleman et al:

    The isothermal boundary conditions give rise to a flow that is strongly influenced by wall-normal gradients of mean density and temperature. These gradients are found to cause an enhanced streamwise coherence of the near-wall streaks, but not to seriously invalidate Morkovin's hypothesis: the magnitude of fluctuations of total temperature and especially pressure are much less than their mean values, and consequently the dominant compressibility effect is that due to mean property variations. The [Van Driest transformation](https://doi.org/10.2514/8.1895) is found to be very successful at both Mach numbers, and when properly scaled, statistics are found to agree well with data from incompressible channel flow results.

## References
1. Coleman et al (1995) https://doi.org/10.1017/S0022112095004587
2. Hamzehloo et al (2021) https://doi.org/10.1002/fld.4879