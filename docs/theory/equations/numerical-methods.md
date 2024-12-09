# Numerical Methods

## HLLC Riemann Solver

Cerisse employs the Harten-Lax-van Leer Contact (HLLC) solver, developed by Toro et al. (1994). The HLLC solver enhances approximate Riemann Solvers by incorporating the intermediate contact wave, a feature that is particularly crucial for accurate modelling in reactive flows applications. To compute the fluxes, the solver is coupled with a second-order Total Variation Diminishing (TVD) scheme, which ensures numerical stability and reduces spurious oscillations

The value of $$U_{RP}$$  in the interface is defined by $$U_l$$  if  $$S_L>0$$,

&#x20; $$U_l^\ast$$ if  $$S_L \leq 0 < S_M$$ , $$U_r^{\ast}$$   if  $$S_M \leq 0  \leq S_R$$   and  $$U_r$$  if $$S_R <0$$

&#x20;The corresponding flux $$F_{RP}$$ is just $$F_{RP}=F(U_{RP})$$.

### The Average State

Assuming the sonic waves speeds, $$S_L$$ and $$S_R$$ are known, we need $$S_M$$ to estimate the average states $$U^\ast$$. The normal velocity and the pressure do not change across a contact discontinuity, therefore the normal velocity is the contact wave speed

$$
S_M=u_l^{\ast}=u_r^{\ast}=u^{\ast}
$$

and the pressure

$$
p_l^{\ast}=p_r^{\ast}=p^{\ast}
$$

where $$q = u$$ is the is the velocity normal to the discontinuity (in this case $u$). The region between the sonic waves has constant pressure $$p^{\ast}$$ and normal velocity $$u^{\ast}$$. To calculate the value of the contact wave speed, the Euler equations across the Riemann fan should be solved, resulting in

$$
S_M = \frac{\rho_r q_r(S_R - q_r)-\rho_l q_l(S_L - q_l) + p_l - p_r}{\rho_r (S_R - q_r)-\rho_l (S_L - q_l)}
$$

From the contact wave speed, we may derive the particle velocity and apply the Rankine-Hugoniot conditions to each acoustic wave to find the average state. In the left wave, the jump relations are

$$
F_l^{\ast} - F_l= S_L (U_l^{\ast} - U_l)
$$

In three dimensional systems the previous equation expand to

$$
S_L \left(\begin{array}{l} \rho_l^{\ast} \\ \rho_l^{\ast}u_l^{\ast} \\ \rho_l^{\ast} v_l^{\ast} \\ \rho_l^{\ast} w_l^{\ast} \\e_l^{\ast} \end{array} \right) - \left(\begin{array}{l} \rho_l^{\ast} q_l^{\ast} \\ \rho_l^{\ast} u_l^{\ast} q_l^{\ast} + p^{\ast}n_x \\ \rho_l^{\ast} v_l^{\ast} q_l^{\ast} + p^{\ast}n_y \\ \rho_l^{\ast} w_l^{\ast} q_l^{\ast} + p^{\ast}n_z \\(e_l^{\ast}+p^{\ast})q_l^{\ast}\end{array} \right) = S_L \left(\begin{array}{l}
\rho_l \\ \rho_l u_l \\ \rho_l v_l\\ \rho_l w_l\\e_l \end{array}\right)- \left( \begin{array}{l} \rho_l q_l \\ \rho_l u_l q_l +p_l n_x\\ \rho_l v_l q_l + p_l n_y \\ \rho_l w_l q_l + p_l n_z \\ (e_l + p_l) q_l
\end{array} \right)
$$

which is easily solved when the contact wave speed is known. Straightforwardly, the first equation gives

$$
\rho_l^{\ast} = \rho_l \frac{S_L-q_l}{S_L-S_M}
$$

Combining both expressions we obtain the intermediate pressure

$$
p^{\ast} = \rho_l (q_L-S_l)(q_l-S_M) + p_l
$$

and solutions to the intermediate-left state (in conserved variables)

$$
(\rho u)_l^{\ast} = \frac{(S_L-q_l)\rho u_l +(p^{\ast}- p_l)}{S_L-S_M}
$$

and

$$
e_l^{\ast} = \frac{(S_L-q_l)\rho E_l -p_lq_l+p^{\ast}S_M}{S_L-S_M}
$$

The procedure to compute the intermediate-right state solution is analogous, but applying the Rankine-Hugoniot condition in the right sonic wave, that is interchanging the subscripts $$l$$ or $$L$$ to $$r$$ and $$R$$, respectively. in previous equations From the two intermediate states the flux may be obtained and replaced in the numerical scheme.

### The Sonic Wave Speed Estimates

To compute the intermediate states the sonic wave speeds ($$S_L,S_R$$) are needed. Following Batten, the wave speeds can be obtained from

$$
S_L = \mbox{min}(q_l - c_l,\tilde{q}-\tilde{c}) \\
S_R = \mbox{max}(q_r + c_r,\tilde{q}+\tilde{c})
$$

where $$\tilde{q}=\tilde{u} n_x+\tilde{v} n_y +\tilde{w} n_z$$ and the average states are given by

$$
\tilde{u}=\frac{(u_l+r_{\rho} u_r)}{(1+r_{\rho})}
$$

The average sound speed is obtained from the average state

$$
\tilde{c}=\sqrt{(\gamma-1)C_p \tilde{T}}
$$

where $$\tilde{T}$$ is obtained through the average enthalpy

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

## Rusanov Scheme

Compact Scheme

## Skew-symmetric

Numerical errors associated with discretisation can be categorised into truncation and aliasing errors (Kravchenko and Moin, 1997[^1]; Lilly, 1965[^2]). Therefore, the concept of numerical order alone is insufficient to fully characterize performance. Key properties such as dissipation, dispersion, and conservation are strongly influenced by the discretization scheme used for the convective term. To illustrate this, consider a one-dimensional scalar equation and three possible formulations for the nonlinear, hyperbolic term:

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
H_{skew} = \frac{1}{2} \left( H_{conv}+ H_{div} \right) 
= \frac{1}{2} \frac{\partial UV }{\partial x} +\frac{1}{2} \left(
U \frac{\partial V }{\partial x} + V \frac{\partial U }{\partial x} \right)
$$

Although the three forms above are equivalent at the continuous level, their discretizations differ significantly in terms of intrinsic properties and performance.

The skew-symmetric form, when used with centered schemes, has been demonstrated to conserve quadratic quantities of interest—such as kinetic energy in the incompressible limit . This conservation is attributed to the reduction of aliasing errors. Furthermore, a Fourier analysis of the three forms reveals that the skew-symmetric formulation possesses superior built-in de-aliasing characteristics (Blaisdell et al. 1996) .

The method implemented here is the approach of Ducros et al. (2000) , which capitalises on the built in de-aliasing property of the skew-symmetric operator of centred schemes, while ensuring local conservation by employing the flux-based formulation:

The flux can be derived in a convective an pressure term

$$
F = UV +F_p = F^{adv}  + F_p
$$

A second order

$$
F^{adv,skew}_{i+1/2} = \frac{1}{4} \left( U_{i} + U_{i+1} \right) \left( V_{i} + V_{i+1}   \right)
$$

and a fourth-order formulation is

$$
F^{adv,skew}_{i+1/2}  =   \frac{1}{3} \left( U_{i} + U_{i+1} \right) \left( V_{i} + V_{i+1}   \right) 
  -\frac{1}{24} \left( U_{i-1} V_{i-1} +   U_{i-1} V_{i+1} +   U_{i} V_{i} +  U_{i} V_{i+2}  +
                            U_{i+1} V_{i-1} +   U_{i+1} V_{i+1} +   U_{i+2} V_{i} +  U_{i+2} V_{i+2}   \right)
$$

At its core, the skew-symmetric scheme, despite its built-in de-aliasing properties, remains a derivative of the family of centered schemes. As such, it can still face stability challenges. Additional strategies to mitigate oscillations and effectively handle shock waves follows

### Artificial dissipation

A way to stabilise the mechanism is through two terms (following Jameson)

$$
\frac{\partial F}{\partial x} \approx \left .\frac{\partial F}{\partial x} \right |_{num}   +
 \frac{\partial^2 \alpha_2 U}{\partial x^2}  + \frac{\partial^4 \alpha_4 U}{\partial x^4}
$$

The second-derivative term is used to capture discontinuities (hereafter _shock_ term) and the fourth-derivative is to control high-frequency noise (_damping_ term). The shock term acts near discontinuities and the damping in smooth parts of the flow. Using the same conservative form as before, the flux is modified by

$$
F^{adv*}_{i+1/2} = F^{skew}_{i+1/2}  +  \left . \alpha_2 \frac{\partial U}{\partial x} \right |_{i+1/2}   +\left . \  \alpha_4 \frac{\partial^3  U}{\partial x^3} \right |_{i+1/2}
$$

The terms can be rewritten using differences, the shock term is :

$$
\left . \alpha_2 \frac{\partial U}{\partial x} \right |_{i+1/2}   \approx \alpha_2 \frac{\Delta U_{i+1/2} }{\Delta x}
$$

where $$\Delta U_{i+1/2} = U_{i+1}-U_i$$. The flux modification is then (grouping constants)

$$
F^{shock}_{i+1/2}  =   \epsilon_{i+1/2}^{(2)} \left( U_{i+1} - U_i \right)
$$

The damping term is similarly

$$
\left .  \alpha_4 \frac{\partial^3  U}{\partial x^3} \right |_{i+1/2}    \approx \frac{\alpha_4}{\Delta x}  \left . \frac{\partial^2  \Delta U}{\partial x^2} \right |_{i+1/2}
$$

Using central differences for $${\partial^2 \Delta U}/{\partial x^2}$$ as

$$
\frac{\partial^2  \Delta U}{\partial x^2} = \frac{\Delta U_{i+3/2} - 2 \Delta U_{i+1/2} + \Delta U_{i-1/2} }{\Delta x^2}   + \mathcal{O}(\Delta x^2)
$$

and replacing the difference, we get

$$
F^{damp}_{i+1/2}  =   \epsilon_{i+1/2}^{(4)}   \left( U_{i+2} - 3 U_{i+1} + 3 U_i - U_{i-1}  \right)
$$

For high-order schemes, the damping term accuracy can be increased, by using a high-order central scheme:

$$
\frac{\partial^2  \Delta U}{\partial x^2} = \frac{ -1 /12 \Delta U_{i+5/2} + 4/3 \Delta U_{i+3/2} - 5/2 \Delta U_{i+1/2} + 4/3 \Delta U_{i-1/2} -1/12 \Delta U_{i-3/2}   }{\Delta x^2}   + \mathcal{O}(\Delta x^4)
$$

the flux is similarly written as

$$
F^{damp}_{i+1/2}  =  \frac{ \epsilon_{i+1/2}^{(4)}  }{12} \left( -U_{i+3} + 17 U_{i+2}  - 46 U_{i+1} +  46 U_i - 17 U_{i-1}  + U_{i-2} \right)
$$

The parameters $$\epsilon_{i+1/2}^{(2)}$$ and $$\epsilon_{i+1/2}^{(4)}$$ control the second and fourth order dissipation.

$$
\epsilon_{i+1/2}^{(2)}=  k^{(2)} | \lambda_{i+1/2} | \psi_{i+1/2}
$$

where $$\psi$$ is the shock/discontinuty detector (1 close to jumps) and $$\lambda$$ is the eigenvalue

$$
\epsilon_{i+1/2}^{(4)}=  \max \left (0, k^{(4)} | \lambda_{i+1/2}|  - \epsilon_{i+1/2}^{(2)}  \right)
$$

The $$- \epsilon_{i+1/2}^{(2)}$$ term, switched off the damping term in the presence of a shock. The shock term is only active if $$\psi > 0$$ The constant values are $$k^{(2)} \sim 0.1-1$$ and $$k^{(4)} \sim 0.01-0.05$$.

### Shock Sensor

The switch $$\psi$$, is computed as the maximum between the sensor value at the node just before and after the interface

$$
\psi_{i+1/2} = \mbox{max}(\psi_i,\psi_{i+1} )
$$

The sensor based on a variable $$\phi$$ is :

$$
\psi_{i} = 2 \frac{|\phi_{i+1} - 2\phi_i + \phi_{i-1} |}{P_{JST} + P_{TVD} +   \varepsilon}
$$

$$\varepsilon$$ is simply a small offset to ensure the denominator is never zero, while $$P_{TVD} = |\phi_{i+1} -\phi_{i} | + |\phi_{i} -\phi_{i-1} |$$and $$P_{JST} = \phi_{i+1} - 2\phi_i + \phi_{i-1}$$

The original formulation works with a sensor on pressure. Cerisse implements the improved approach of Bouheraoua (2014) by including an additional density sensor and coupling the two as :

$$
\psi = \frac{\psi_\rho^2 + \psi_P^2}{\psi_\rho + \psi_P}
$$

## WENO and TENO

Weighted Essentially Non-Oscillatory (WENO) methods, introduced by Liu et al. (1994), employ a nonlinear adaptive procedure to automatically select the locally smoothest stencil. This approach aims to avoid using stencils that cross discontinuities when interpolating the interface flux.

The WENO family encompasses various variations, which can be further classified. Despite these differences, all WENO methods share a common feature: the interface flux is expressed as a linear combination of fluxes derived from the stencils.

$$
F_{i+1/2} = \sum_s w_s F^s_{i+2}  \;\; \;\;  w_s = \frac{\alpha_s}{\sum_s \alpha_s}
$$

Reconstruction in characteristic variables improves performance, as the post-shock oscillations are reduced. WENO is known to be excessively dissipative in smooth parts of the flow

### TENO

Designed[^3] to reduce numerical dissipation further than WENO

## KEEP

Central KEEP (_non-dissipative and physically-consistent kinetic energy and entropy preserving_) schemes for compressible flows[^4]

### References

[^1]: Kravchenko, A. and Moin, P. (1997). On the effect of numerical errors in large eddy simulations of turbulent flows. [Journal of Computational physics, 131(2):310–322](https://doi.org/10.1006/jcph.1996.5597)

[^2]: Lilly, D. K. (1965). On the computational stability of numerical solutions of time-dependent non-linear geophysical fluid dynamics problems. [Monthly Weather Review, 93(1):11–25](https://doi.org/10.1175/1520-0493\(1965\)093%3C0011:OTCSON%3E2.3.CO;2)

[^3]: Fu, L., Hu, X. Y., and Adams, N. A. (2017). Targeted eno schemes with tailored resolution property for hyperbolic conservation laws. [Journal of Computational Physics, 349:97–121](https://doi.org/10.1016/j.jcp.2017.07.054)

[^4]: Yuichi Kuya, Soshi Kawai, (2020) A stable and non-dissipative kinetic energy and entropy preserving (KEEP) scheme for non-conforming block boundaries on Cartesian grids, [Computers & Fluids, Volume 200, 104427](https://doi.org/10.1016/j.compfluid.2020.104427)
