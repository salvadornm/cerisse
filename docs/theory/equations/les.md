---
icon: filter
---

# LES

## Filtered Navier-Stokes Equations

Filter process The spatial filter of a function $$\phi$$ is defined as its convolution integral with a filter function $$G$$, according to:

$$
\overline{\phi} = \int_V  \phi(\mathbf{x}',t)  G(\mathbf{x}-\mathbf{x}',\Delta) dV'
$$

where $$\Delta$$ is the characteristic filter width in each respective direction and $$\overline{\phi}$$ is the filtered quantity.

The actual form of the filter is usually _implicit_ and is not needed while solving the LES equations. Is common to use a _box/top-hat filter_ of the form below, where the cut-off scale/filter width is taken proportional to the mesh size $$\Delta \propto h$$

$$
G(\mathbf{x}-\mathbf{x}',\Delta) =
\left \{ 
\begin{array}{lc}
{1}/{\Delta^3} & |\mathbf{x}-\mathbf{x}'| < \Delta/2\\
0 & \mbox{otherwise}\\
\end{array}
\right.
$$

Other filter kernels are possible (see [Pope's book](https://www.cambridge.org/highereducation/books/turbulent-flows/C58EFF59AF9B81AE6CFAC9ED16486B3A#overview) for more detail on LES filters). In the finite volume method, when the **filter width matches the local cell size** , _i.e._, $$\Delta= h$$ and the **filter used is the box filte**r; the _cell-averaged_ value of a variable is equivalent to its _filtered_ value.

For variable density flows it is convenient to introduce the mass-weighted Favre filtering operation :

$$
\widetilde{\phi } = \frac{\overline{\rho \phi} }{\bar{\rho}}
$$

It is common to assume _commutability_ of the filtering and derivative operators, that is to say

$$
\overline{ \frac{\partial \phi }{\partial x_j} } = \frac{\partial \overline{\phi}}{\partial x_j}
$$

This relation is only true under several assumptions, including the restrictive condition of constant filter width throughout the domain. Is common to neglect this error and assume the effects are incorporated into the sub-grid model.

#### Continuity

Applying the filtering operator to the [continuity equation](dns.md)

$$
\frac{\partial \overline{\rho} }{\partial t} + \frac{\partial \overline{\rho u_j} }{\partial x_j} = 0
$$

and using Favre weigthing

$$
\frac{\partial \bar{\rho} }{\partial t} + \frac{\partial \bar{\rho} \tilde{u}_j }{\partial x_j} = 0
$$

#### Filtered Momentum

$$
\frac{\partial \overline{\rho u_i} }{\partial t} + \frac{\partial 
    \overline{\rho u_i u_j}}{\partial x_j} = -\frac{\partial \bar{p}}{\partial x_i} + \frac{\partial \overline{\tau}_{ij}}{\partial x_j}
$$

Using Favre weigthing

$$
\frac{\partial \bar{\rho} \widetilde{u}_i }{\partial t} + \frac{\partial \bar{\rho} \widetilde{u_i u_j}}{\partial x_j} = -\frac{\partial \bar{p}}{\partial x_i} + \frac{\partial \overline{\tau}_{ij}}{\partial x_j}
$$

#### Filtered Species transport

$$
\frac{\partial \overline{\rho_k} }{\partial t} + \frac{\partial 
\overline{ \rho_k u_j }}{\partial x_j} = \frac{\partial \overline{J}_j^k}{\partial x_j} +  \overline{\rho \dot{\omega}_k}
$$

Using Favre weigthing

$$
\frac{\partial \bar{\rho} \bar{Y}_k }{\partial t} + \frac{\partial \bar{\rho}
\widetilde{ u_j Y_k}}{\partial x_j} = \frac{\partial \overline{J}_j^k}{\partial x_j} + \bar{\rho} \widetilde{\dot{\omega}_k}
$$

#### Filtered Energy

$$
\frac{\partial  \overline{E_t}  }{\partial t} + \frac{\partial \overline{(E_t + p) u_j}}{\partial x_j} = \frac{\overline{\partial u_i \tau_{ij}}}{\partial x_j} - 
\frac{\partial \overline{q}_j }{\partial x_j} + \sum_k \frac{ \overline{\partial h_k J_j^k}}{\partial x_j}
$$

## Sub-grid Closures

The filtered momentum equation requires modelling of the following terms

$$
\overline{\rho u_i u_j} =\bar{\rho} \widetilde{u_i u_j}  = \bar{\rho} \tilde{u}_i \tilde{u}_j  - \tau^{sgs}_{ij}
$$

where a sub-gris stress tensor is introduced that need modelling

$$
\tau^{sgs}_{ij} =   \bar{\rho} \tilde{u}_i \tilde{u}_j  - \bar{\rho} \widetilde{u_i u_j}
$$

A common strategy found in many LES studies is the use of an eddy viscosity-like assumption founded on _Boussinesq’s hypthothesis_

$$
\tau_{ij}^{sgs} - \frac{1}{3} \tau_{kk} = \mu_{sgs} 
\left( \tilde{S}_{ij} - \frac{1}{3} \tilde{S}_{kk} \right)
$$

A sub-grid viscosity is introduced, similar to turbulent viscosity in RANS-type models that requires modelling. Implemented models include:

#### Smagorinsky Model

The Smagorinsky model [\[3\]](les.md#references) assumes that small, unresolved turbulent eddies behave like an **eddy viscosity**, enhancing momentum diffusion, with a characteristic length scale $$\ell_{sgs}$$.&#x20;

$$
\mu_{sgs}  = \bar{\rho} (C_S \Delta) ^2 || \tilde{S}_{ij} ||
$$

with $$\tilde{S}_{ij}$$represents the filtered strain tensor and $$C_S$$ the **Smagorinsky constant,** with values between **0.1** and **0.2**.\
&#x20;$$||\tilde{S}_{ij}|| = \sqrt{2 \tilde{S}_{ij} \tilde{S}_{ij} }$$ is the Frobenius norm of the filtered strain tensor. The length scale $$l_{sgs}= C_S \Delta$$ is a sub-grid length scale, which can be consider proportional to the integral length-scale $$\ell$$. The Smagorinsky model (and most gradeint models) assumes that small sub-grid scales are isotropic.

#### WALE Model

The eddy viscosity in the WALE model [\[1\]](les.md#references) is computed

$$
\mu_{sgs}  = \bar{\rho} (C_w \Delta) ^2  
 \frac{(\mathcal{S}_{ij}^d\mathcal{S}_{ij}^d)^{3/2}}{(S_{ij} S_{ij})^{5/2} - (\mathcal{S}_{ij}^d \mathcal{S}_{ij}^d)^{5/4}}
$$

where $$\mathcal{S}_{ij}^d = \mathcal{S}_{ij} - 1/3 \mathcal{S}_{kk}$$ is the traceless, symmetric tensor of the square of the velocity gradient

$$
\mathcal{S}_{ij} = \frac{1}{2} \left(\frac{\partial u_i}{\partial x_k} \frac{\partial u_k}{\partial x_j} + \frac{\partial u_j}{\partial x_k} \frac{\partial u_k}{\partial x_i} \right)
$$

Model coefficient is in the range $$C_w = \sqrt{10.6} \cdot C_S = 0.325 - 0.5$$

#### Diffusivity and heat

The sub-grid transport of a scalar is splitted in

$$
\widetilde{u_j Y_k} =  \tilde{u}_j \tilde{Y}_k   -  D_{sgs} \frac{\partial \tilde{Y}_k}{\partial x_j}
$$

where $$D_{sgs}$$ is the sub-grid diffusivity, which taken proportional to sub-grid viscosity

$$
\bar{\rho}  D_{sgs} = \frac{\mu_{sgs}}{\text{Sc}_{sgs}}
$$

with $$\text{Sc}_{sgs}$$ is a constant often taken  between **0.4** and **1**, or found using a dynamic procedure. All species diffuse at the smallest scales _due to turbulence_ at the same speed.

$$
\overline{\rho u_j (e_t + P/\rho )}  = \bar{\rho} \widetilde{ u_j h_t} = \bar{\rho}\tilde{u}_j \tilde{h}_t - \lambda_{sgs} \frac{\partial \tilde{T}}{\partial x_j}
$$

where $$h_t \equiv  e_t + P/\rho$$is the specific total enthalpy and $$\lambda_{sgs}$$ a sub-grid conductivity, that can be related to the sub-grid viscosity through:

$$
\lambda_{sgs} = \frac{\mu_{sgs} C_p}{\text{Pr}_{sgs}}
$$

Where the sub-grid Prandtl number is introduced, which is a constant taken in the range **0.4-1**. Cerisse works with the ratio to molecular Prandtl number and redefines the  sub-grid conductivity as

$$
\frac{\lambda_{sgs}}{\lambda} = \frac{\mu_{sgs}}{\mu} \frac{\text{Pr}}{\text{Pr}_{sgs}}
$$

To avoid computing the specific heat.  The ratio $$\text{Pr}/\text{Pr}_{sgs}$$ is often less than 1 in gases, with a common choice of 0.5/0.7  $$\approx$$**0.7**  . In general, if $${\text{Pr}}/{\text{Pr}_{sgs}} > 1$$, sub-grid turbulent eddies transport momentum more efficiently than heat (_vs. molecular case_). If $${\text{Pr}}/{\text{Pr}_{sgs}} < 1$$, sub-grid turbulent eddies transport momentum more efficiently than heat (_vs. molecular case_).\
\
The isotropic part of the sub-grid stress $$⅓ \tau_{kk}$$ is neglected in incompressible flows (absorved in the pressure)  and is often modelled using **Yoshizawa** model [\[2\]](les.md#references) in compressible flows&#x20;

$$
\tau_{kk} = \bar{\rho} C_I \Delta^2 ||\tilde{S}_{ij}||^2
$$

where  $$C_I$$ is a model constant taken often as  **0.008**. This expression can be used to estimate the sub-grid kinetic energy

$$
k_{sgs} =\frac{3}{2} \tau_{kk}  = C_Y \Delta^2 ||\tilde{S}_{ij}||^2
$$

with $$C_Y$$ taken as **0.0066** .

### Other LES  stategies

Other LES approaches do not use gradient-type models. For example, **Implicit Large Eddy Simulation** (ILES) is a class of turbulence modelling techniques where **no explicit sub-grid scale  model is added**. Instead, the **numerical discretisation itself acts as the turbulence model**. ILES assumes that a carefully chosen numerical scheme (e.g., monotonicity-preserving, shock-capturing, or flux-limited) can do this naturally. ILES is very sensitive to mesh and numerics,  and it may show more numerical dissipation at small scales[\[4\]](les.md#references). However, with  appropriate schemes,  they can recover the correct total energy decay rate. Similarly, using explicit sub-grid models close to shocks is not recommended.\
To select this model in Cerisse  ran DNS with appropiate numerics (for example [WENO-type schemes](numerical-methods.md)).\
Alternatively, joint velocity LES-PDF models [\[5\]](les.md#references) do not need gradient-type sub-grid models (see [Turbulent Combustion modelling](turbcomb.md) for a more detailed description)

### Other unknowns

In conventional LES, fluctuations of transport properties are assumed to be small within the filter width and therefore, molecular fluxes can be approximated by

$$
\overline{q}_j \approx - \lambda (\tilde{T}) \frac{\partial \tilde{T}}{\partial x_j}
$$

Molecular fluxes scale with the inverse of Reynolds number, $$\text{Re}^{-1}$$, making them relatively small in turbulent flows. Consequently, errors associated with molecular transport properties often (but not always) have a small impact in the solution of turbulent flows. &#x20;

#### References

\[1] Nicoud, F., Ducros, F. (1999) Subgrid-Scale Stress Modelling Based on the Square of the Velocity Gradient Tensor. [_Flow, Turbulence and Combustion_ 62, 183–200](https://doi.org/10.1023/A:1009995426001)\
\[2] Yoshizawa, A.   Horiuti,  K. A. (1985) Statistically-Derived Subgrid-Scale Kinetic Energy Model for the Large-Eddy Simulation of Turbulent Flows.[ _Journal of the Physical Society of Japan, 54, 2834-2839_](https://doi.org/10.1143/JPSJ.54.2834)\
\[3] Smagorinsky, J. (1963). General circulation experiments with the primitive equations. [_Monthly Weather Review, 91(3):99–164._](https://doi.org/10.1175/1520-0493\(1963\)091%3C0099:GCEWTP%3E2.3.CO;2)\
\[4] Garnier, E , Mossi, M , Sagaut , P, Comte, P, Deville, M (1999) On the Use of Shock-Capturing Schemes for Large-Eddy Simulation. [_Journal of Computational Physics 153(2): 273-311_](https://doi.org/10.1006/jcph.1999.6268)\
\[5] Un, T-H, Navarro-Martinez, S. Stochastic fields with adaptive mesh refinement for high-speed turbulent combustion, [_Combustion and  Flame_, 272, 113897 (2025)](https://doi.org/10.1016/j.combustflame.2024.113897) \
\[6] Boris J. P, Grinstein, F. F. , Oran, E. S. Kolbe, R. S. (1992) New insights into large eddy simulation. [_Fluid Dynamics Research 10(4): 199-228_](https://doi.org/10.1016/0169-5983\(92\)90023-P)



\
\
