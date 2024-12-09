---
icon: filter
---

# LES

## Filtered Navier-Stokes Equations

Filter process The spatial filter of a function $$\phi$$ is defined as its convolution integral with a filter function $$G$$, according to:

$$
\overline{\phi} = \int_V  \phi(\mathbf{x}',t)  G(\mathbf{x}-\mathbf{x}',\Delta) dV'
$$

where $$\Delta$$ is the characteristic filter width in each respective direction and $$\overline{\phi}$$is the filtered quantity.

The actual form of the filter is usually _implicit_ and not needed while soling the LES equations. Is common to use a box/top-hat filter of the form below where the cut-off scale/filter width is usually taken proportional not the mesh size $$\Delta \propto h$$

$$
G(\mathbf{x}-\mathbf{x}',\Delta) =
\left \{ 
\begin{array}{lc}
{1}/{\Delta^3} & |\mathbf{x}-\mathbf{x}'| < \Delta/2\\
0 & \mbox{otherwise}\\
\end{array}
\right.
$$

Where other filter kernels are possible (see Pope's book). In the finite volume method, when the filter width matches the local cell size , i.e., $$\Delta= h$$and the filter used is the box filter (as described above), the cell-averaged value of a variable is equivalent to its filtered value.

For variable density flows it is convenient to introduce the mass-weighted Favre filtering operation :

$$
\widetilde{\phi } = \frac{\overline{\rho \phi} }{\bar{\rho}}
$$

It is commonb to assume assumption of commutability of the filtering and derivative operators, that is to say

$$
\overline{ \frac{\partial \phi }{\partial x_j} } = \frac{\partial \overline{\phi}}{\partial x_j}
$$

This relation is only true under several assumptions, including the restrictive condition of constant filter width throughout the domain. Is common to neglect this error and assume the effects are incorporated into the SGS models.

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

The filtered momentum equation requires modelling for the SGS stress tensor.

$$
\tau_{ij}^{SGS}= \overline{\rho u_i u_j} =\bar{\rho} \widetilde{u_i u_j}
$$

A common strategy found in many LES studies is the use of an eddy viscosity-like assumption founded on Boussinesq’s hypthothesis

$$
\tau_{ij}^{sgs} - \frac{1}{3} \tau_{kk} = \mu_{sgs} 
\left( \tilde{S}_{ij} - \frac{1}{3} \tilde{S}_{kk} \right)
$$

Where a sub-grid viscosity is introduced.

#### Smagorinsky

$$
\mu_{sgs}  = \bar{\rho} (C_S \Delta) ^2 || \tilde{S}_{ij} ||
$$

with $$\tilde{S}_{ij}$$represents the filtered strain tensor and $$C_S$$the Smagorinsky constant, with values between 0.1-0.2. $$||\tilde{S}_{ij}|| = \sqrt{2 \tilde{S}{ij} \tilde{S}_{ij} }$$ is the Frobenius norm of the filtered strain tensor. The length scale $$l_{sgs}= C_S \Delta$$ is a sub-grid lengths scale, which can be consider proportional to the integral lenght-scale $$\ell$$.

#### WALE

The eddy viscosity in the WALE model is computed

$$
\mu_{sgs}  = \bar{\rho} (C_w \Delta) ^2  
 \frac{(\mathcal{S}_{ij}^d\mathcal{S}_{ij}^d)^{3/2}}{(S_{ij} S_{ij})^{5/2} - (\mathcal{S}_{ij}^d \mathcal{S}_{ij}^d)^{5/4}}
$$

where $$\mathcal{S}_{ij}^d = \mathcal{S}_{ij} - 1/3 \mathcal{S}_{kk}$$ is the traceless, symmetric tensor of the square of the velocity gradient

$$
\mathcal{S}_{ij} = \frac{1}{2} \left(\frac{\partial u_i}{\partial x_k} \frac{\partial u_k}{\partial x_j} + \frac{\partial u_j}{\partial x_k} \frac{\partial u_k}{\partial x_i} \right)
$$

Model coefficients are in the range $$C_w = 0.325 - 0.5$$

#### Diffusivty and heat

The sub-grid transport of a scalar is splitted in

$$
\widetilde{u_j Y_k} =  \tilde{u}_j \tilde{Y}_k   -  D_{sgs} \frac{\partial \tilde{Y}_k}{\partial x_j}
$$

where $$D_{sgs}$$ is the sub-grid diffusivity, which taken proportional to sub-grid viscosity

$$
\bar{\rho}  D_{sgs} = \frac{\mu_{sgs}}{Sc_{sgs}}
$$

with $$Sc_{sgs}$$ a constant often taken as 0.4-1.

In conventional LES, fluctuations of transport properties are assumed to be small within the filter width and therefore, molecular fluxes can be approximated by

$$
\overline{q}_j \approx - \lambda (\tilde{T}) \frac{\partial \tilde{T}}{\partial x_j}
$$

Molecular fluxes scale with $$Re^{-1}$$, making them relatively small in turbulent flows. Consequently, errors associated with molecular transport properties often have a minimal impact.

#### References
