# LES Documentation


## Governing Equations

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


### Mass Fraction of specie k

$$
\frac{\partial \rho Y_k }{\partial t} + \frac{\partial \rho u_j Y_k}{\partial x_j}   
=   \frac{\partial J_j^k}{\partial x_j} +  \rho \dot{\omega}_k
$$


### Total Energy

$$
\frac{\partial E_t }{\partial t} + \frac{\partial (E_t + p) u_j}{\partial x_j}   = 
 \frac{\partial u_i \tau_{ij}}{\partial x_j}+
\frac{\partial q_j }{\partial x_j}
+  \sum_k  \frac{ \partial h_k J_j^k}{\partial x_j}
$$


## Closures


### Transport fluxes

Newtonian flows, Fourier and Fickian diffusion and also include (optionally) Soret diffusion


$$
q_j = -\lambda  \frac{\partial T }{\partial x_j}
$$

$$
\tau_{ij} = 2 \mu S_{ij}
$$

$$
J_j^k = - D_k \frac{\partial X_k }{\partial x_j}
$$


**Mixing Rules**



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

Assiming the filter function conmutes with spatial derivatives


### Continuity

$$
\frac{\partial \overline{\rho}  }{\partial t}   
+ \frac{\partial \overline{\rho u} }{\partial x}
+ \frac{\partial \overline{ \rho v} }{\partial y}
+ \frac{\partial \rho w }{\partial z}
= 0
$$


or using Favre (or density weighted) filtering
$$
\overline{\rho \phi} = \bar{\rho} \widetilde{\phi} 
$$

we obtain

$$
\frac{\partial \bar{\rho}  }{\partial t}   
+ \frac{\partial \bar{\rho} \widetilde{u} }{\partial x}
+ \frac{\partial \bar{ \rho} \widetilde{v} }{\partial y}
+ \frac{\partial \bar{\rho} \widetilde{ w} }{\partial z}
= 0
$$

### Momentum

$$
\frac{\partial \overline{\rho u_i} }{\partial t} + \frac{\partial \overline{\rho u_i u_j}}{\partial x_j}   = -\frac{\partial \bar{p}}{\partial x_i} + \frac{\partial \overline{\tau}_{ij} }{\partial x_j}
$$

or


$$
\frac{\partial \bar{\rho} \widetilde{u_i} }{\partial t} + \frac{\partial \bar{\rho} \widetilde{ u_i u_j}}{\partial x_j}   = -\frac{\partial \bar{p}}{\partial x_i} + \frac{\partial \overline{\tau}_{ij} }{\partial x_j}
$$

## Conventional Sub-grid Closures



