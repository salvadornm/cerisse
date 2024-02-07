# PDF Documentation

LES-PDF mode


## Stochastic Fields Equations


The final equations to be solved are

### Continuity

$$
 \frac{d \varrho^n}{d t} + \frac{\partial \varrho^n \mathscr{U}^n_j }{\partial x_j} = 0
 \label{finalcontinuitySPDE}
$$


### Momentum

$$
\frac{d \varrho^n \mathscr{U}^n_i}{d t} + \frac{\partial \varrho^n \mathscr{U}^n_j \mathscr{U}^n_i}{\partial x_j} = -\frac{\partial \mathscr{P}^n}{\partial x_i} + \frac{\varrho^n}{\overline{\rho}} \frac{\partial \widetilde{\tau}_{ij}}{\partial x_i} + \varrho^n G_{ij} \left( \mathscr{U}^n_j - \widetilde{u}_j \right) + \varrho^n \sqrt{C_0 \frac{\epsilon_{sgs}}{\overline{\rho}}} \frac{d W^n_i}{d t}
$$

### Mass Fraction of specie k

$$
\frac{\partial \varrho^n \mathscr{Y}^n_k}{\partial t} +\frac{\partial  \varrho^n \mathscr{U}^n_i \mathscr{Y}^n_k}{\partial x_i} =  \frac{\varrho^n}{\overline{\rho}}\frac{\partial \widetilde{J}_{k,i} }{\partial x_i} + \varrho^n \dot{\omega}_k - \frac{1}{2} C_{Y_k} \frac{\epsilon_{sgs}}{k_{sgs}} \varrho^n \left( \mathscr{Y}^n_k - \widetilde{Y}_k \right)
$$ 

###  Energy
$$
\frac{\partial \varrho^n \mathscr{E}^n_t}{\partial t} + \frac{\partial \varrho^n \mathscr{U}^n_i \mathscr{E}^n_t}{\partial x_i} = \frac{\varrho^n}{\overline{\rho}}\frac{\partial \widetilde{q}_i}{\partial x_i}  - \frac{\varrho^n}{\overline{\rho}}\frac{\partial \overline{p} \widetilde{u}_i}{\partial x_i} + \frac{\varrho^n}{\overline{\rho}}\frac{\partial \widetilde{\tau}_{ij}\widetilde{u}_j}{\partial x_i}  - \frac{1}{2} C_{e_t} \frac{\epsilon_{sgs}}{k_{sgs}} \varrho^n \left( \mathscr{E}^n_t - \widetilde{e_t} \right)
$$


The employed closure relation for the dissipation of the sub-grid kinetic energy is

$$
\epsilon_{sgs} = C_\epsilon k_{sgs}^{3/2}/\Delta
$$

where the constant $$C_\epsilon = 1.05 $$. 
The micro-mixing constants $$C_{Y_k} = 2 $$ 
 and the Langevin constant is set to 2.1,


$$
k_{sgs}= \frac{1}{2} \left( \frac{1}{N_f} \sum_{n=1}^{N_f} (\mathscr{U}^n_i-\widetilde{u}_i)^2 \right)
$$

The filtered variables can be obtained from the average of the Eulerian stochastic fields.

$$
\overline{Q}  = \frac{1}{N_f}\sum_{n=1}^{N_f} Q^n; \ \ \ \ \widetilde{Q} = \frac{\sum_{n=1}^{N_f} \varrho^n Q^n}{\sum_{n=1}^{N_f} \varrho^n}
$$



## Implementation Notes





