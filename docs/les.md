# LES Documentation


## Governing Equations

### Continuity

$$
\frac{\partial \rho}{\partial t} + \frac{\partial \rho u_j}{\partial x_j} = 0
$$

### Momentum

$$
\frac{\partial \rho u_i}{\partial t} + \frac{\partial \rho u_i u_j}{\partial x_j} = -\frac{\partial p}{\partial x_i} + \frac{\partial \tau_{ij}}{\partial x_j}
$$

### Mass Fraction of species $k$

$$
\frac{\partial \rho Y_k}{\partial t} + \frac{\partial \rho u_j Y_k}{\partial x_j} = \frac{\partial J_j^k}{\partial x_j} +  \rho \dot{\omega}_k
$$

### Total Energy

$$
\frac{\partial E}{\partial t} + \frac{\partial (E + p) u_j}{\partial x_j} = \frac{\partial u_i \tau_{ij}}{\partial x_j} + \frac{\partial q_j}{\partial x_j}
$$


## Closures

### Transport fluxes

Newtonian flows, Fourier heat flux, and Hirschfelder-Curtiss
approximation diffusion are used:

$$
q_j = -\lambda \frac{\partial T }{\partial x_j} + \sum_{k} h_k J_j^k
$$

$$
\tau_{ij} = 2 \mu S_{ij} + \left(\mu_b - \frac{2}{3}\mu\right) S_{kk} \delta_{ij}
$$

$$
J_j^k = \rho D_k \frac{Y_k}{X_k} \left(\frac{\partial X_k}{\partial x_i} + (X_k-Y_k)\frac{\partial \ln p}{\partial x_i}\right)
$$


### Equations of state

$$
p = f(\rho,T,Y_k)
$$

Perfect gas, ideal gas, Van der Waals and Soave-Redlich-Kwong

The calorific equation of state is

$$
e_i =  f(\rho,T,Y_k)
$$


## Filtered Equations

Assuming that the filter function conmutes with spatial derivatives,

### Continuity

$$
\frac{\partial \overline{\rho}  }{\partial t}   
+ \frac{\partial \overline{\rho u_i} }{\partial x_i}
= 0
$$
or using Favre (density-weighted) filtering $\overline{\rho \phi} = \bar{\rho} \widetilde{\phi}$, we obtain

$$
\frac{\partial \bar{\rho}}{\partial t}   
+ \frac{\partial \bar{\rho} \widetilde{u}_i}{\partial x_i}
= 0
$$

### Momentum

$$
\frac{\partial \bar{\rho} \widetilde{u}_i}{\partial t} + \frac{\partial \bar{\rho} \widetilde{u}_i \widetilde{u}_j}{\partial x_j} = -\frac{\partial \bar{p}}{\partial x_i} + \frac{\partial \overline{\tau}_{ij}}{\partial x_j} - \frac{\partial \tau^{sgs}_{ij}}{\partial x_j}
$$

### Species

$$
\frac{\partial \bar{\rho} \widetilde{Y}_k }{\partial t} + \frac{\partial \bar{\rho} \widetilde{u}_j \widetilde{Y}_k}{\partial x_j}   
= \frac{\partial \bar{J}_j^k}{\partial x_j} + \frac{\partial J^{sgs}_{k, j}}{\partial x_j} +  \overline{\rho \dot{\omega}_k}
$$

### Energy

$$
\frac{\partial \bar{\rho} \widetilde{E}}{\partial t} + \frac{\partial \bar{\rho} \widetilde{E} \widetilde{u}_j}{\partial x_j} = -\frac{\partial \bar{p} \widetilde{u}_i}{\partial x_i} + \frac{\partial \overline{\tau}_{ij} \widetilde{u}_i}{\partial x_j} + \frac{\partial \bar{q}_{j}}{\partial x_j} + \frac{\partial q^{sgs}_{j}}{\partial x_j}
$$


## Conventional Sub-grid Closures

Using the Boussinesq hypothesis,

$$
  \tau^{sgs}_{ij} \approx -2 \bar{\rho} \nu_{sgs} \left( \tilde{S}_{ij} - \frac{\delta_{ij}}{3} \tilde{S}_{kk} \right)
$$
$$
  \frac{\partial J^{sgs}_{k, j}}{\partial x_j} \approx \frac{\partial}{\partial x_j} \left( -\bar\rho \frac{\nu_{sgs}}{{Sc}_{sgs}} \frac{\partial \tilde Y_k}{\partial x_j} \right)
$$
$$
  \frac{\partial q^{sgs}_j}{\partial x_j} \approx \frac{\partial}{\partial x_j} \left( -\bar\rho c_p \frac{\nu_{sgs}}{{Pr}_{sgs}} \frac{\partial T}{\partial x_j} \right)
$$

### Smagorinsky model
$$
  \nu_{sgs} = (C_s \Delta)^2 |\tilde{S}| \,, \qquad    
  |\tilde{S}| = (2\,\tilde{S}_{ij}\tilde{S}_{ij})^{1/2}
$$
where $C_s$ is the Smagorinsky constant, usually taken between 0.1 and 0.2.

### WALE model
The wall-adapting local eddy-viscosity (WALE) model is designed to yield the correct asymptotic scaling of $\nu_{sgs} \propto y^3$ near walls, where $y$ is wall normal distance.
$$
  \nu_{sgs} = (C_w \Delta)^2 \frac{(S_{ij}^d S_{ij}^d)^{3/2}}{(\tilde S_{ij} \tilde S_{ij})^{3/2} + (S_{ij}^d S_{ij}^d)^{5/2}}
$$
$$
  S_{ij}^d = \frac{1}{2}\left(\frac{\partial\tilde u_i}{\partial x_k}\frac{\partial\tilde u_k}{\partial x_j} + \frac{\partial\tilde u_j}{\partial x_k}\frac{\partial\tilde u_k}{\partial x_i} \right) + \frac{1}{3} \frac{\partial\tilde u_l}{\partial x_k}\frac{\partial\tilde u_k}{\partial x_l} \delta_{ij}
$$    
with $C_w \approx 0.5$. Compared to the Smagorinsky model, the WALE model ensures correct wall shear stress prediction and produces zero $\nu_{sgs}$ in pure two-dimensional shear flows.


