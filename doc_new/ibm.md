# Inmesred Boundaries


## Theory


RHS adds..

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
