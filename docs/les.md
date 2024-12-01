# Filtered Navier-Stokes Equations

Filter process
The spatial filter of a function $$\phi$$ 
is defined as its convolution integral with a filter
function $$G$$, according to:

$$
\overline{\phi} = \int_V  \phi(\mathbf{x}',t)  G(\mathbf{x}-\mathbf{x}',\Delta) dV'
$$

where $$\Delta$$ is the characteristic filter width in
each respective direction and $$\overline{\phi} $$ is the filtered quantity.

The actual form of the filter is usually *implicit* and not need. Is common 
to use a box/top-hat filter of the form below where the cut-off scale/filter
width is usually taken as proprtional not the mesh size $$\Delta \propto h$$

$$
G(\mathbf{x}-\mathbf{x}',\Delta) =
\left \{ 
\begin{array}{lc}
{1}/{\Delta^3} & |\mathbf{x}-\mathbf{x}'| < \Delta/2\\
0 & \mbox{otherwise}\\
\end{array}
\right.
$$

Where other filter kernels are possible (see Pope's book).
In the finite volume method, when the filter width 
matches the local cell size , i.e., $$\Delta= h $$
and the filter used is the box filter (as described above), the cell-averaged value of a variable is equivalent to its filtered value.

For variable density flows it is convenient to introduce the mass-weighted Favre filtering
operation :
$$
\widetilde{\phi } = \frac{\overline{\rho \phi} }{\bar{\rho}}
$$

It is commonb to assume assumption of commutability of
the filtering and derivative operators, that is to say

$$
\overline{ \frac{\partial \phi }{\partial x_j} } = \frac{\partial \overline{\phi}}{\partial x_j} 
$$

This relation is only true under several assumptions, including the  restrictive condition
of constant filter width throughout the domain.
Is common to  neglect this error and assume the effects are incorporated
into the SGS models.

### Continuity

Applying the filtering operator to the continuty equation (link to DNS)

$$
\frac{\partial \overline{\rho} }{\partial t} + \frac{\partial \overline{\rho u_j} }{\partial x_j} = 0
$$

and using Favre weigthing

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






### References

Omer Rathore, "Numerical Simulation of Combustion
Instability" *PhD Thesis*, Imperial College London (2022)
