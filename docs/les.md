# Filtered Navier-Stokes Equations

Filter process

$$
\overline{\phi} = \int_V G(r-r',\Delta) dV'
$$

Favre weighting

$$
\widetilde{\phi } = \frac{\overline{\rho \phi} }{\bar{\rho}}
$$

### Continuity

$$
\frac{\partial \overline{\rho} }{\partial t} + \frac{\partial \overline{\rho u_j} }{\partial x_j} = 0
$$

Using Favre weigthing

$$
\frac{\partial \bar{\rho} }{\partial t} + \frac{\partial \bar{\rho} \tilde{u}_j }{\partial x_j} = 0
$$

### Filtered Momentum

$$
    \frac{\partial \overline{\rho u_i} }{\partial t} + \frac{\partial 
    \overline{\rho u_i u_j}}{\partial x_j} = -\frac{\partial \bar{p}}{\partial x_i} + \frac{\partial \overline{\tau}_{ij}}{\partial x_j}
$$


Using Favre weigthing
$$
 \frac{\partial \bar{\rho} \widetilde{u}_i }{\partial t} + \frac{\partial \bar{\rho} \widetilde{u_i u_j}}{\partial x_j} = -\frac{\partial \bar{p}}{\partial x_i} + \frac{\partial \overline{\tau}_{ij}}{\partial x_j}
 $$

### Filtered Species transport

$$
\frac{\partial \overline{\rho_k} }{\partial t} + \frac{\partial 
\overline{ \rho_k u_j }}{\partial x_j} = \frac{\partial \overline{J}_j^k}{\partial x_j} +  \overline{\rho \dot{\omega}_k}
$$

Using Favre weigthing
$$
\frac{\partial \bar{\rho} \bar{Y}_k }{\partial t} + \frac{\partial \bar{\rho}
\widetilde{ u_j Y_k}}{\partial x_j} = \frac{\partial \overline{J}_j^k}{\partial x_j} + \bar{\rho} \widetilde{\dot{\omega}_k}
$$

### Filtered Energy


