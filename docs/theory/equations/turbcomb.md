---
icon: fire
---

# Turbulent Combustion

Turbulent Combustion model

## ATF

The Artificial Thickening Method (ATF) approach essentially involves "thickening" the flame, allowing the mesh to resolve the scalar gradients across it. This serves a dual purpose: it enables capturing the flame dynamics and reduces numerical diffusion errors associated with resolving the sharp scalar jumps in the flame. The aim of ATF methods is to achieve adequate resolution within the flame front (_i.e._, the region with the sharpest gradients) on the transformed mesh, even when the coarse LES mesh lacks sufficient resolution. The geometry transformation ratio, denoted as $$\mathcal{F}$$, is referred to as the thickening factor.

$$
\frac{\partial \rho Y_k}{\partial \tau} +\frac{\partial \rho u_i Y_k}{\partial \xi_i}= \frac{\partial}{\partial \xi_i }\left(\rho D \mathcal{F}\frac{\partial Y_k}{\partial \xi_i}\right) + \frac{ \dot{\omega}}{\mathcal{F}}
$$

The transformation "thickens" the flame but maintains the correct flame speed by increasing diffusion accordingly.

The sub-grid wrinkling requires modelling and this is usually done in the ATF context by the inclusion of a so called efficiency function. The efficiency function, $$E$$, is defined by a dimensionless wrinkling factor $$\mathcal{E}$$, and its ratio between a laminar flame compared with its thickened counterpart. Inclusion into the ATF model results in a modification to the reactive scalar transport equation:

$$
\frac{\partial \rho Y_k}{\partial \tau} +\frac{\partial \rho u_i Y_k}{\partial \xi_i}= \frac{\partial}{\partial \xi_i }\left(\rho D E\mathcal{F}\frac{\partial Y_k}{\partial \xi_i}\right) + \frac{ E\dot{\omega}}{\mathcal{F}}
$$

$$
E=\frac{\mathcal{E}\vert_{\delta_f=\delta_f^0}}{\mathcal{E}\vert_{\delta_f=\hat{\delta}_f^0}} \geq 1
$$

The factor $$\mathcal{E}$$ essentially relates the the total flame front wrinkling with its resolved component. It can be approximated as:

$$
\mathcal{E} \approx 1+\beta \Delta \vert <\nabla  \vec{n} >_{sgs}\vert
$$

where $$\beta$$ is a constant and $$\vert <\nabla \vec{n} >_{sgs}\vert$$ is the SGS surface curvature.

## PasR

The Partially Stired Reactor Model

## Eulerian Stochastic Fields

The Stochastic fields equations for the joint-velocity-scalar energy PDF equations to be solved are

### Continuity

$$
\frac{d \varrho^n}{d t} + \frac{\partial \varrho^n \mathscr{U}^n_j }{\partial x_j} = 0
$$

### Momentum

$$
\frac{d \varrho^n \mathscr{U}^n_i}{d t} + \frac{\partial \varrho^n \mathscr{U}^n_j \mathscr{U}^n_i}{\partial x_j} = -\frac{\partial \mathscr{P}^n}{\partial x_i} + \frac{\varrho^n}{\overline{\rho}} \frac{\partial \widetilde{\tau}_{ij}}{\partial x_i} + \varrho^n G_{ij} \left( \mathscr{U}^n_j - \widetilde{u}_j \right) + \varrho^n \sqrt{C_0 \frac{\epsilon_{sgs}}{\overline{\rho}}} \frac{d W^n_i}{d t}
$$

### Mass Fraction of specie k

$$
\frac{\partial \varrho^n \mathscr{Y}^n_k}{\partial t} +\frac{\partial  \varrho^n \mathscr{U}^n_i \mathscr{Y}^n_k}{\partial x_i} =  \frac{\varrho^n}{\overline{\rho}}\frac{\partial \widetilde{J}_{k,i} }{\partial x_i} + \varrho^n \dot{\omega}_k - \frac{1}{2} C_{Y} \frac{\epsilon_{sgs}}{k_{sgs}} \varrho^n \left( \mathscr{Y}^n_k - \widetilde{Y}_k \right)
$$

### Energy

$$
\frac{\partial \varrho^n \mathscr{E}^n_t}{\partial t} + \frac{\partial \varrho^n \mathscr{U}^n_i \mathscr{E}^n_t}{\partial x_i} = \frac{\varrho^n}{\overline{\rho}}\frac{\partial \widetilde{q}_i}{\partial x_i}  - \frac{\varrho^n}{\overline{\rho}}\frac{\partial \overline{p} \widetilde{u}_i}{\partial x_i} + \frac{\varrho^n}{\overline{\rho}}\frac{\partial \widetilde{\tau}_{ij}\widetilde{u}_j}{\partial x_i}  - \frac{1}{2} C_{E} \frac{\epsilon_{sgs}}{k_{sgs}} \varrho^n \left( \mathscr{E}^n_t - \widetilde{e_t} \right)
$$

The employed closure relation for the dissipation of the sub-grid kinetic energy is

$$
\epsilon_{sgs} = C_\epsilon k_{sgs}^{3/2}/\Delta
$$

where the constant $$C_\epsilon = 1.05$$. The micro-mixing constants $$C_Y = 2$$ and the Langevin constant is set to 2.1,

$$
k_{sgs}= \frac{1}{2} \left( \frac{1}{N_f} \sum_{n=1}^{N_f} (\mathscr{U}^n_i-\widetilde{u}_i)^2 \right)
$$

The filtered variables can be obtained from the average of the Eulerian stochastic fields.

$$
\overline{\phi}  = \frac{1}{N_f}\sum_{n=1}^{N_f} \phi^n; \ \ \ \ \widetilde{\phi} = \frac{\sum_{n=1}^{N_f} \varrho^n \phi^n}{\sum_{n=1}^{N_f} \varrho^n}
$$

### References

1. Omer Rathore, "Numerical Simulation of Combustion Instability" _PhD Thesis_, Imperial College London (2022)&#x20;
2. Yuri Almeida "Large Eddy Simulation of Supersonic Combustion using a Probability Density Function method" _PhD Thesis_, Imperial College London (2019)
3. Tin-Hang Un and Salvador Navarro-Martinez, “Stochastic fields with adaptive mesh refinement for high-speed turbulent combustion”, [_Comb. Flame_, 272, 113897 (2025)](https://doi.org/10.1016/j.combustflame.2024.113897)
