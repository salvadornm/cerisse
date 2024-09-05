# Immersed Boundaries


## Theory

In *Immersed Boundary Methods * (IBM), boundary conditions are enforced by introducing source terms,  using interpolation  to impose the influence of the boundary on the fluid flow.

## Set-up input


Define STL file with geometry

```
#-------------------------------- IMMERSED BOUNDARY ---------------------------#
ib.filename = sphere_fine_20k.stl
```


## Set-up prob.h


```cpp
pp.add   ("ib.move",0); // 0=false, 1=true
pp.add   ("ib.plot_surf",0); // 0=false, 1=true
```

## Set-up bcnormal



## Implementation: IBM class

```eib_t``` is explicit geometry (triangulation based) immersed boundary method
class. It holds an array of IBMultiFab, one for each AMR level; and it also holds
the geometry
