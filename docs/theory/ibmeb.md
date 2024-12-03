# Geometry

Embedded/Immersed boundary methods enable the representation of complex geometries on structured grids. Unlike body-fitted approaches, which require intricate coordinate transformations, these methods involve minimal modifications to the flux derivative calculations. Additionally, they are better suited for fluid-structure interaction problems compared to body-fitted methods.
Embedded boundary methods are generally categorized into two main types: cut-cell methods (EBM) and ghost-point methods (IBM). Similar to other numerical techniques, each approach has its own advantages and disadvantages.


## Immersed Boundaries

Immersed Boundaries, or Ghost-Point Immersed Boundary Methods (GPIBM), originate from finite difference discretization and are not inherently designed to ensure discrete conservation at the embedded boundary.
These methods approximate the nearest point to the embedded boundary, known as the ghost point/cell, in a way that ensures the correct boundary flux is satisfied.

FIGURE

## Embedded Boundaries


Embedded Boundaries, or Cut-Cell Methods, originate from finite volume discretization and are specifically designed to ensure discrete conservation at the embedded boundary.
These methods represent the boundary using a piecewise reconstruction within a Cartesian mesh, resulting in a sharp and accurate depiction of the interface. Instead of approximating a cell value, they focus on approximating the surface flux.

FIGURE


### References

[1]: Monal Patel, "Hypersonic Flows Around Complex
Geometries with Adaptive Mesh Refinement and Immersed Boundary Method" *PhD Thesis*, Imperial College London (2022)
