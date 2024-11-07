# Equations

## Navier-Stokes Equations

### Continuity

$$
\frac{\partial \rho  }{\partial t}   
+ \frac{\partial \rho u_j }{\partial x_j}
= 0
$$

### Momentum

$$
\frac{\partial \rho u_i }{\partial t} + \frac{\partial \rho u_i u_j}{\partial x_j}   = -\frac{\partial p}{\partial x_i} + \frac{\partial \tau_{ij}}{\partial x_j}
$$


### Chemical species transport k

$$
\frac{\partial \rho Y_k }{\partial t} + \frac{\partial \rho u_j Y_k}{\partial x_j}   
=   \frac{\partial J_j^k}{\partial x_j} +  \rho \dot{\omega}_k
$$

The code solves for the partial densities $$\rho_k = \rho Y_k$$ 

$$
\frac{\rho_k }{\partial t} + \frac{ \rho_k u_j }{\partial x_j}   
=   \frac{\partial J_j^k}{\partial x_j} +  \rho \dot{\omega}_k
$$


### Total Energy

$$
\frac{\partial E_t }{\partial t} + \frac{\partial (E_t + p) u_j}{\partial x_j}   = 
 \frac{\partial u_i \tau_{ij}}{\partial x_j}+
\frac{\partial q_j }{\partial x_j}
+  \sum_k  \frac{ \partial h_k J_j^k}{\partial x_j}
$$

where $$E_t = \rho e_i + \rho k $$ is the total energy plus kinetic energy.

## Closures


### Transport fluxes

The main cosures assume Newtonian flows, Fourier heat transfer and and Fickian diffusion and also  (optionally) Soret diffusion
The **Shear stress** is given by

$$
\tau_{ij} = 2 \mu S_{ij}  + (\mu_b - \frac{2}{3} \mu)\nabla \cdot \vec{v} \delta_{ij}
$$

where $$\mu$$ and $$\mu_b$$ are the viscosity and bulk viscosity respectively.
The **Heat Flux**

$$
q_j = -\lambda  \frac{\partial T }{\partial x_j}
$$

where $$\lambda$$ is the conductivity.
And the **diffusion**

$$
J_j^k = - D_k \frac{\partial X_k }{\partial x_j}
$$


**Mixing Rules**



### Equations of state

The equations of state are of the form

$$
p = f(\rho,T,Y_k)
$$

and the calorific equation of state

$$
e_i =  f(\rho,T,Y_k)
$$

The current implementations include  perfect gas, ideal gas and Soave-Redlich-Kwong (through PelePhysics).





