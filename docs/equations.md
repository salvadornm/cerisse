# Equations

Cerisse solves equations in the general form

$$
\frac{\partial \mathbf{U} }{\partial t} + \frac{\partial \mathbf{F}_j }{\partial x_j}  + \frac{\partial \mathbf{F}_{visc} }{\partial x_j} =  \mathbf{S}
$$

Where $$\mathbf{U}$$ is the array in conservative variables

$$
\mathbf{U} = \begin{bmatrix} \rho \\ \rho u_i  \\ E_t \\ \rho Y_1 \\ \rho Y_2 \\ ... \end{bmatrix}
$$

and $$\mathbf{F}$$ represent the fluxes of the conserved variables across a surface.
For example in x-direction.

$$
\mathbf{F}_x =
\begin{bmatrix}
\rho u \\ \rho u_2 + p \\ \rho u v \\ \rho u w \\ u(E_t + p) \\ \rho u Y_1 \\ \rho u Y_2 \\ ...
\end{bmatrix}
$$
These are called the The Euler fluxes. If these fluzes are the only present, the resultant equations are the **Euler Equations** 

$$\mathbf{F}_{visc}$$ represents the fluxes due to molecular transport (viscous streess, heat fluxes, mass disffusion, etc).
While $$\mathbf{S}$$ is a generic source term.
While combining these terms appropiately we can build differnt type of equations (check PROB), such as Euler, Navier-Stokes, Reactive Navier-Stokes, Non-equilbrium thermodynamics, etc.
The reactive equations can be see in the DNS tab.

## Finite Volume Method

The expression can be written in divergence form

$$
\frac{\partial \mathbf{U} }{\partial t} =  
-\nabla \cdot \left( \mathbf{F} +  \mathbf{F}_{visc} \right) + \mathbf{S} =  -\nabla \cdot \mathbf{F}^\ast + \mathbf{S}
$$

Integrating over cell and dividing over control volume $V_{ijk}$
$$
\frac{1}{V_{ijk}} \int \frac{\partial \mathbf{U} }{\partial t} dV =  
\frac{1}{V_{ijk}}  \int -\nabla \cdot \mathbf{F}^\ast  \, dV +  \frac{1}{V_{ijk}}   \int \mathbf{S} dV
$$

Defining the cell-averaged value as
$$
\phi_{ijk} = \frac{1}{V_{ijk}}  \int \phi dV   
$$
and using divergence theorem
$$
\frac{\partial \mathbf{U}_{ijk} }{\partial t} =  
\frac{1}{V_{ijk}}  \oint_{\delta V} \mathbf{F}^\ast  \cdot \mathbf{n} \, dA +  \mathbf{S}_{ijk}
$$
Using polyhedral cells, the following expression  follows by summing over the cells faces
$$
 \frac{\partial \mathbf{U}_{ijk} }{\partial t}  =  
\frac{1}{V_{ijk}}  \sum_{f} \mathbf{F}_f^\ast \cdot \mathbf{n}_f \, A_f +  \mathbf{S}_{ijk}
$$
The above formulation ensures global conservation as the fluxes of
 shared faces between two control volumes cancel out.