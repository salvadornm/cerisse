# Geometry

Embedded/Immersed boundary methods enable the representation of complex geometries on structured grids. Unlike body-fitted approaches, which require intricate coordinate transformations, these methods involve minimal modifications to the flux derivative calculations. Additionally, they are better suited for fluid-structure interaction problems compared to body-fitted methods.
Embedded boundary methods are generally categorized into two main types: cut-cell methods (EBM) and ghost-point methods (IBM). Similar to other numerical techniques, each approach has its own advantages and disadvantages.


## Immersed Boundaries

Immersed Boundaries, or Ghost-Point Immersed Boundary Methods (GPIBM), originate from finite difference discretization and are not inherently designed to ensure discrete conservation at the embedded boundary.
These methods approximate the nearest point to the embedded boundary, known as the ghost point/cell, in a way that ensures the correct boundary flux is satisfied.

FIGURE IBM

## Embedded Boundaries


Embedded Boundaries, or Cut-Cell Methods, originate from finite volume discretization and are specifically designed to ensure discrete conservation at the embedded boundary.
These methods represent the boundary using a piecewise reconstruction within a Cartesian mesh, resulting in a sharp and accurate depiction of the interface. Instead of approximating a cell value, they focus on approximating the surface flux.

FIGURE EBM


In EB the flux divergence at cell next to the wall (cut-cell) is

$$
\nabla \cdot F = \sum_{k=1}^{Nfaces} \frac{F_k A _k}{V} = 
\frac{F_n A_n - F_s A_s +  F_e A_e - F A_w + F_r A_r -F_l A_l }{ V} + F_{wall}{V}
$$

The volume of the cells is $$V= \phi h^3$$ (assuming isotropic mesh), where $$\phi$$  is the "filled" fraction of the cell.
Using that the area of the face is proportional to the mesh size  
$$A_k = f_k h^2$$. Where $$f_k$$ is the ratio of obunsutcred area with the maximum. It will be 1 if the area is unobstructed.

Reaplacing

$$
\nabla \cdot F = \sum_{k=1}^{6} \frac{F_k f _k}{ \phi h} + \frac{F_{wall}}{\phi h}
$$

where the flux at the wall would be defined based on BC type.

### References

[1]: Monal Patel, "Hypersonic Flows Around Complex
Geometries with Adaptive Mesh Refinement and Immersed Boundary Method" *PhD Thesis*, Imperial College London (2022)

