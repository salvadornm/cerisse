# Numerical Methods

Cerisse has implemented several numerical methods within the finite volume framework. While most of these methods are designed for high-speed flows, not all are limited to such cases. The method discussed here specifically addresses the Euler (convective) term of the equations. In compressible flows, high-frequency noise can accumulate when physical scales are not adequately resolved, especially near discontinuities. To mitigate this, many methods either directly filter or dissipate the noise, or use upwind-type stencils that indirectly introduce dissipation. It is important to note that no perfect numerical scheme exists; a compromise must be made between accuracy, speed, and stability. In the sections below, a brief overview of the rationale behind the numerical methods is provided. This description is not exhaustive, and the reader is encouraged to refer to the original papers for a more comprehensive explanation.

## HLLC Riemann Solver

Cerisse employs the Harten-Lax-van Leer Contact (HLLC) solver, developed by Toro et al. (1994). The HLLC solver enhances approximate Riemann Solvers by incorporating the intermediate contact wave, a feature that is particularly crucial for accurate modelling in reactive flows applications.

<figure><img src="../../.gitbook/assets/rieman.jpg" alt=""><figcaption><p>Scheme of the HLLC Riemann Solver with three waves propagating ar speed SL, SM and SR.<br>Two acoustic waves and one contact. The waves separate four constant states <span class="math">U_L</span> , <span class="math">U_L^\ast</span>, <span class="math">U_R^\ast</span>, <span class="math"> U_R</span></p></figcaption></figure>

&#x20;The value of $$U_{RP}$$ in the interface is defined by

$$
U_{RP}= \left\{ \begin{array}{cl}
        U_L & \mbox{if $S_L >0$} \\
	U_L^{\ast} & \mbox{if $S_L \leq 0 < S_M$}\\
	U_R^{\ast} & \mbox{if $S_M \leq 0  \leq S_R$} \\
	U_R & \mbox{if $S_R < 0$}
           \end{array}
\right.
$$

The corresponding flux $$F_{RP}$$ is just $$F_{RP}=F(U_{RP})$$.

### The Average State

Assuming the sonic waves speeds, $$S_L$$ and $$S_R$$ are known, we need $$S_M$$ to estimate the intermediate average states $$U^\ast$$. The normal velocity and the pressure do not change across a contact discontinuity (mechanical equilbrium) and therefore the normal velocity is the contact wave speed

$$
S_M=u_l^{\ast}=u_r^{\ast}=u^{\ast}
$$

and the pressure

$$
p_l^{\ast}=p_r^{\ast}=p^{\ast}
$$

where $$q = u$$ is the is the velocity _normal_ to the discontinuity. The region between the sonic waves has constant pressure $$p^{\ast}$$ and normal velocity $$u^{\ast}$$. To calculate the value of the contact wave speed, the Euler equations across the Riemann fan should be solved, resulting in

$$
S_M = \frac{\rho_r q_r(S_R - q_r)-\rho_l q_l(S_L - q_l) + p_l - p_r}{\rho_r (S_R - q_r)-\rho_l (S_L - q_l)}
$$

From the contact wave speed,  the Rankine-Hugoniot conditions are applied to each acoustic wave to find the average state.  In the left wave, the jump relations are:

$$
F_l^\ast - F_l = S_L (U_l^\ast - U_l)
$$

The density of the inytermediate left state is:&#x20;

$$
\rho_l^\ast = \rho_l \frac{S_L - q_l}{S_L - S_M }
$$

and solutions to the intermediate-left state (in conserved variables)

$$
(\rho u)_l^{\ast} = \frac{(S_L-q_l)\rho u_l +(p^{\ast}- p_l)}{S_L-S_M}
$$

and

$$
e_l^{\ast} = \frac{(S_L-q_l)\rho E_l -p_lq_l+p^{\ast}S_M}{S_L-S_M}
$$

The procedure to compute the intermediate-right state solution is analogous, but applying the Rankine-Hugoniot condition in the right sonic wave, that is interchanging the subscripts $$l$$ or $$L$$ to $$r$$ and $$R$$, respectively. in previous equations. From the two intermediate states the flux may be obtained and replaced in the numerical scheme. The first-order scheme is

$$
F_{i+1/2}= F(U_{RP})
$$

### The Sonic Wave Speed Estimates

To compute the intermediate states the sonic wave speeds ($$S_L,S_R$$) are needed. Following Batten, the wave speeds can be obtained from

$$
S_L = \mbox{min}(q_l - c_l,\tilde{q}-\tilde{c})
$$

and

$$
S_R = \mbox{min}(q_r + c_r,\tilde{q} + \tilde{c})
$$

where $$\tilde{q}=\tilde{u} n_x+\tilde{v} n_y +\tilde{w} n_z$$ and the average state are (same for otehr components)



where $$\tilde{T}$$ is obtained through the average enthalpy

$$
\tilde{u}=\frac{(u_l+r_{\rho} u_r)}{(1+r_{\rho})}
$$

$$
\tilde{H}=\frac{(H_l+r_{\rho} H_r)}{(1+r_{\rho})}
$$

and

$$
\tilde{H}=C_p \tilde{T} + \frac{1}{2} \left( \tilde{u}^2+\tilde{v}^2+\tilde{w}^2 \right)
$$

where $$r_{\rho}$$ is the ratio of densities

$$
r_{\rho}=\sqrt{\rho_r/\rho_l}
$$

### High order extension

The above scheme is formally first-order, a high order extension can be build by using a Total Variation Diminishing (TVD) reconstruction

## Rusanov Scheme

The Rusanov flux is a simple upwind flux that requires a single wave-speed estimate. In the current implementation, it is very compact and can be used to perform quick tests. The numerical flux function is typically given by:

$$
F_{i+1/2}= \frac{1}{2} \left( F_{i+1} + F_{i} \right) - \frac{1}{2} \left( \left| \lambda_{i+1} \right| + \left| \lambda_i \right| \right) (U_{i+1} - U_i)
$$



where $$\lambda_i$$ ​ is the local characteristic speed at the i-th cell (often taken as the maximum eigenvalue of the Jacobian of the flux function).

### Skew-symmetric

Numerical errors associated with discretisation can be categorised into truncation and aliasing errors (Kravchenko and Moin, 1997; Lilly, 1965). Therefore, the concept of numerical order alone is insufficient to fully characterize performance. Key properties such as dissipation, dispersion, and conservation are strongly influenced by the discretization scheme used for the convective term. To illustrate this, consider a one-dimensional scalar equation and three possible formulations for the nonlinear, hyperbolic term:

$$
\frac{\partial U}{\partial t} + \frac{\partial H U}{\partial x} = 0
$$

$$
H_{div} =\frac{\partial UV }{\partial x}
$$

$$
H_{conv} = U \frac{\partial V }{\partial x} + V \frac{\partial U }{\partial x}
$$

$$
H_{skew} = \frac{1}{2} \left( H_{conv}+ H_{div} \right) = \frac{1}{2} \frac{\partial UV }{\partial x} +\frac{1}{2} \left( U \frac{\partial V }{\partial x} + V \frac{\partial U }{\partial x} \right)
$$



Although the three forms above are equivalent at the continuous level, their discretizations differ significantly in terms of intrinsic properties and performance.

The skew-symmetric form, when used with centered schemes, has been demonstrated to conserve quadratic quantities of interest—such as kinetic energy in the incompressible limit . This conservation is attributed to the reduction of aliasing errors. Furthermore, a Fourier analysis of the three forms reveals that the skew-symmetric formulation possesses superior built-in de-aliasing characteristics (Blaisdell et al. 1996) .

The method implemented here is the approach of Ducros et al. (2000) , which capitalises on the built in de-aliasing property of the skew-symmetric operator of centred schemes, while ensuring local conservation by employing the flux-based formulation:

The flux can be derived in a convective an pressure term

$$
F = UV +F_p = F^{adv} + F_p
$$

A second order scheme can be then

$$
F^{adv,skew}_{i+1/2} = \frac{1}{4} \left( U{i} + U_{i+1} \right) \left( V_{i} + V_{i+1} \right)
$$

and a fourth-order formulation scheme, despite its built-in de-aliasing properties, remains a derivative of the family of centered schemes. As such, it can still face stability challenges. Additional strategies to mitigate oscillations and effectively handle shock waves follows

#### Artificial dissipation

A way to stabilise the mechanism is through two terms (following Jameson)

$$
\frac{\partial F}{\partial x} \approx \left .\frac{\partial F}{\partial x} \right |_{num} + \frac{\partial^2 \alpha_2 U}{\partial x^2} + \frac{\partial^4 \alpha_4 U}{\partial x^4}
$$

The second-derivative term is used to capture discontinuities (hereafter _shock_ term) and the fourth-derivative is to control high-frequency noise (_damping_ term). The shock term acts near discontinuities and the damping in smooth parts of the flow. Using the same conservative form as before, the flux is modified by

$$
F^{adv}_{i+1/2} = F^{skew}_{i+1/2} + \left . \alpha_2 \frac{\partial U}{\partial x} \right |_{i+1/2} +\left . \ \alpha_4 \frac{\partial^3 U}{\partial x^3} \right |{i+1/2}
$$

The terms can be rewritten using differences, the shock term is&#x20;

$$
\left . \alpha_2 \frac{\partial U}{\partial x} \right |_{i+1/2} \approx \alpha_2 \frac{\Delta U_{i+1/2} }{\Delta x}
$$

where $$\Delta U_{i+1/2} = U_{i+1}-U_i$$. The flux modification is then (grouping constants)

$$
F^{shock}_{i+1/2} = \epsilon_{i+1/2}^{(2)} \left( U_{i+1} - U_i \right)
$$

The damping term is similarly

&#x20;$$\left . \alpha_4 \frac{\partial^3 U}{\partial x^3} \right |_{i+1/2} \approx \frac{\alpha_4}{\Delta x} \left . \frac{\partial^2 \Delta U}{\partial x^2} \right |_{i+1/2}$$

Using central differences for $${\partial^2 \Delta U}/{\partial x^2}$$ as

$$
\frac{\partial^2 \Delta U}{\partial x^2} = \frac{\Delta U_{i+3/2} - 2 \Delta U_{i+1/2} + \Delta U_{i-1/2} }{\Delta x^2} + \mathcal{O}(\Delta x^2)
$$

and replacing the difference, we get

F^{damp}_{i+1/2} = \epsilon_{i+1/2}^{(4)} \left( U\_{i+2} - 3 U\_{i+1} + 3 U\_i - U\_{i-1} \right)

$$
For high-order schemes, the damping term accuracy can be increased, by using a high-order central scheme:
$$

\frac{\partial^2 \Delta U}{\partial x^2} = \frac{ -1 /12 \Delta U\_{i+5/2} + 4/3 \Delta U\_{i+3/2} - 5/2 \Delta U\_{i+1/2} + 4/3 \Delta U\_{i-1/2} -1/12 \Delta U\_{i-3/2} }{\Delta x^2} + \mathcal{O}(\Delta x^4)

$$
the flux is similarly written as
$$

F^{damp}_{i+1/2} = \frac{ \epsilon_{i+1/2}^{(4)} }{12} \left( -U\_{i+3} + 17 U\_{i+2} - 46 U\_{i+1} + 46 U\_i - 17 U\_{i-1} + U\_{i-2} \right)

$$$
The parameters $$\epsilon_{i+1/2}^{(2)}$$ and $$\epsilon_{i+1/2}^{(4)}$$ control the second and fourth order dissipation.
$$$

\epsilon\_{i+1/2}^{(2)}= k^{(2)} | \lambda\_{i+1/2} | \psi\_{i+1/2}

$$$
where $$\psi$$ is the shock/discontinuty detector (1 close to jumps) and $$\lambda$$ is the eigenvalue
$$$

\frac{\partial \rho k }{\partial t} + \frac{\partial \rho k }{\partial x\_j} = - u \frac{\partial p }{\partial x\_j}

$$$
The sensor based on a variable $$\phi$$ is :
$$$

\psi\_{i} = 2 \frac{|\phi\_{i+1} - 2\phi\_i + \phi\_{i-1} |}{P\_{JST} + P\_{TVD} + \varepsilon}

$$$
$$\varepsilon$$ is simply a small offset to ensure the denominator is never zero, while $$P_{TVD} = |\phi_{i+1} -\phi_{i} | + |\phi_{i} -\phi_{i-1} |$$and $$P_{JST} = \phi_{i+1} - 2\phi_i + \phi_{i-1}$$

The original formulation works with a sensor on pressure. Cerisse implements the improved approach of Bouheraoua (2014) by including an additional density sensor and coupling the two as :
$$$

\psi = \frac{\psi\_\rho^2 + \psi\_P^2}{\psi\_\rho + \psi\_P}

$$
## WENO and TENO

Weighted Essentially Non-Oscillatory (WENO) methods, introduced by Liu et al. (1994), employ a nonlinear adaptive procedure to automatically select the locally smoothest stencil. This approach aims to avoid using stencils that cross discontinuities when interpolating the interface flux.

The WENO family encompasses various variations, which can be further classified. Despite these differences, all WENO methods share a common feature: the interface flux is expressed as a linear combination of fluxes derived from the stencils.
$$

F\_{i+1/2} = \sum\_s w\_s F^s\_{i+2} ;; ;; w\_s = \frac{\alpha\_s}{\sum\_s \alpha\_s}

$$
Reconstruction in characteristic variables improves performance, as the post-shock oscillations are reduced. WENO is known to be excessively dissipative in smooth parts of the flow

### TENO

Designed to reduce numerical dissipation further than WENO

## KEEP

Central KEEP (_non-dissipative and physically-consistent kinetic energy and entropy preserving_) schemes for compressible flows
These scheme are base din splitting the energy equation.
$$

\frac{\partial E\_t }{\partial t} + \frac{\partial (\rho e + \rho k + p) u\_j}{\partial x\_j} = 0

$$$
where $$E_t = \rho e + \rho k$$ is the total energy plus kinetic energy.
In the inviscid limit, from the momentum equation is possible to derive the kinetic energy equation
$$$

\frac{\partial \rho k }{\partial t} + \frac{\partial \rho u\_j k }{\partial x\_j} + u \frac{\partial p }{\partial x\_j} = 0

$$
which implies
$$

\frac{\partial \rho e }{\partial t} + \frac{\partial \rho u\_j e }{\partial x\_j} + p \frac{\partial u }{\partial x\_j} = 0

$$
from
$$

d e = T d s - \frac{p}{\rho^2} d \rho

$$
which implies conservation of entropy
$$

\frac{\partial \rho s }{\partial t} + \frac{\partial \rho u\_j s }{\partial x\_j} = 0

$$$
KEEP stil need a shock capturing term.

### Split terms

Given a function $$ f = a b $$, there are two forms to split the derivate

Divergence
$$$

\frac{\partial f}{\partial x} = \frac{\partial a b}{\partial x}

$$
Quadratic
$$

\frac{\partial f}{\partial x} =\
\frac{1}{2}\left( a \frac{\partial b }{\partial x} + b \frac{\partial a }{\partial x} + \frac{\partial a b}{\partial x} \right)

$$$
These formulations are analytically equivalent. With a function $$ f = a b c $$,  there are three forms to split the derivate

Divergence
$$$

\frac{\partial f}{\partial x} = \frac{\partial a b c}{\partial x}

$$
Quadratic
$$

\frac{\partial f}{\partial x} =\
\frac{1}{2}\left( a \frac{\partial b c}{\partial x} + b c \frac{\partial a }{\partial x} + \frac{\partial abc }{\partial x} \right)

$$
Cubic
$$

\frac{\partial f}{\partial x} =\
\frac{1}{4}\left( a \frac{\partial b c}{\partial x} + b \frac{\partial a c }{\partial x} + c \frac{\partial a b}{\partial x} + a b \frac{\partial c }{\partial x} + a c \frac{\partial b }{\partial x} + b c \frac{\partial a }{\partial x} +\
\frac{\partial abc }{\partial x} \right)

$$
### References

[1]: Morinishi, Y. (1995). Conservative properties of finite difference schemes for incompressible flow. 
[Center for Turbulence Research Annual Research Briefs](https://ntrs.nasa.gov/citations/19960022304)

[2]: Blaisdell, G., Spyropoulos, E., and Qin, J. (1996). The effect of the formulation of nonlinear terms on aliasing errors in spectral methods. 
[Applied Numerical Mathematics, 21(3):207–219](https://doi.org/10.1016/0168-9274(96)00005-0)

[3]: Ducros, F., Laporte, F., Soulères, T., Guinot, V., Moinat, P., and Caruelle, B. (2000). High order fluxes for conservative skew-symmetric-like schemes in structured meshes: application to compressible flows. 
[Journal of Computational Physics, 161(1):114–139](https://doi.org/10.1006/jcph.2000.6492)

[4]: Kravchenko, A. and Moin, P. (1997). On the effect of numerical errors in large eddy
simulations of turbulent flows. 
[Journal of Computational physics, 131(2):310–322](https://doi.org/10.1006/jcph.1996.5597)

[5]: Lilly, D. K. (1965). On the computational stability of numerical solutions of time-dependent non-linear geophysical fluid dynamics problems. 
[Monthly Weather Review, 93(1):11–25](https://doi.org/10.1175/1520-0493(1965)093%3C0011:OTCSON%3E2.3.CO;2) 

[6]: Liu, X.-D., Osher, S., and Chan, T. (1994). Weighted essentially non-oscillatory schemes.
[Journal of Computational physics, 115(1):200–212](https://doi.org/10.1006/jcph.1994.1187)

[7]: Yuichi Kuya, Soshi Kawai, (2020) A stable and non-dissipative kinetic energy and entropy preserving (KEEP) scheme for non-conforming block boundaries on Cartesian grids,
[Computers & Fluids, Volume 200, 104427](https://doi.org/10.1016/j.compfluid.2020.104427)

[8]: Fu, L., Hu, X. Y., and Adams, N. A. (2017). Targeted eno schemes with tailored resolution property for hyperbolic conservation laws. 
[Journal of Computational Physics, 349:97–121](https://doi.org/10.1016/j.jcp.2017.07.054)

[9]: Toro, E. F., Spruce, M., and Speares, W. (1994). Restoration of the contact surface in the hll-riemann solver. 
[Shock waves, 4(1):25–34](https://doi.org/10.1007/BF01414629)

[10]: Bouheraoua, L. (2014). Simulation aux grandes échelles et modélisation de la combustion
supersonique. [PhD thesis](https://theses.hal.science/tel-01197487v1), Rouen, INSA.
$$
