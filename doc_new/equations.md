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


### Mass Fraction of specie k



### Total Energy

$$
\frac{\partial E_t }{\partial t} + \frac{\partial (E_t + p) u_j}{\partial x_j}   = 
 \frac{\partial u_i \tau_{ij}}{\partial x_j}+
\frac{\partial q_j }{\partial x_j}
+  \sum_k  \frac{ \partial h_k J_j^k}{\partial x_j}
$$


## Closures


### Transport fluxes

Newtonian flows, Fourier and Fickian diffusion and also including (optionally) Soret diffusion


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

They are of the form

$$
p = f(\rho,T,Y_k)
$$

Current implementations throigh Pele-Physics include 
perfect gas, ideal gas, Van der Waals and Soave-Redlich-Kwong

The calorific equation of state is in the form

$$
e_i =  f(\rho,T,Y_k)
$$



## Filtered Equations for LES

Assiming the filter function conmutes with spatial derivatives


### Filtered Continuity

$$
\frac{\partial \overline{\rho}  }{\partial t}   
+ \frac{\partial \overline{\rho u_j} }{\partial x_j}
= 0
$$


or using Favre (or density weighted) filtering
$$
\overline{\rho \phi} = \overline{\rho} \widetilde{\phi} 
$$

we obtain

$$
\frac{\partial \bar{\rho}  }{\partial t}   
+ \frac{\partial \bar{\rho} \widetilde{u}_j }{\partial x_j}
= 0
$$

### Filtered Momentum

$$
\frac{\partial \overline{\rho u_i} }{\partial t} + \frac{\partial \overline{\rho u_i u_j}}{\partial x_j}   = -\frac{\partial \bar{p}}{\partial x_i} + \frac{\partial \overline{\tau}_{ij} }{\partial x_j}
$$

or


$$
\frac{\partial \bar{\rho} \widetilde{u_i} }{\partial t} + \frac{\partial \bar{\rho} \widetilde{ u_i u_j}}{\partial x_j}   = -\frac{\partial \bar{p}}{\partial x_i} + \frac{\partial \overline{\tau}_{ij} }{\partial x_j}
$$



### Filtered Species Transport


$$
\frac{\partial \bar{\rho} \widetilde{Y}_k }{\partial t} + \frac{\partial \bar{\rho} \widetilde{u_j Y_k}}{\partial x_j}   
=   \frac{\partial \overline{J}_j^k}{\partial x_j} +  \bar{\rho} \widetilde{\dot{\omega}_k}
$$

### Filtered Energy


## Conventional Sub-grid Closures



