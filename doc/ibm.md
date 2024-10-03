## Boundary Condtions

### input

Define in input file apply to all faces in that direction

```
# -1 = N/A (Incase of periodic)
# 0 = Interior           3 = Symmetry
# 1 = Inflow             4 = SlipWall
# 2 = Outflow            5 = NoSlipWall
# 6 = user defined
cns.lo_bc = 1 1
cns.hi_bc = 2 1
```

### bcnormal

During run-time, complex bounday conditions can be set0-up by defining the
```bcnormal``` function in ```prob.h```

```cpp
bcnormal(const Real x[AMREX_SPACEDIM], Real dratio, const Real s_int[5],
         const Real s_refl[ProbClosures::NCONS], Real s_ext[5], const int idir,
         const int sgn, const Real time, GeometryData const &,
         ProbClosures const &closures, ProbParm const &prob_parm)
```

This function is called for every face in the boundary.
The function needs to define the array ```s_ext[5]`` the value at the ghost points. The arguments are


| Option                      | Type          | Dimensions | Description                                                  |
| --------------------------- | ------------- |:-------:| ------------------------------------------------------------ |
| ```x```       | Real          | SPACEDIM   | ghost cell cooordinates  |
| ```dratio```  | Real          | 1          |  ghost/first internal distance ratio                          |
| ```s_int```   | Real          | NCONS      | flow state inside domain  |
| ```s_ext```   | Real          | NCONS      | flow state to be filled  |
| ```idir```    | int           | 1          | direction (0: x, 1: y, 2: z)  |
| ```sign```    | int           | 1          | high or low boundary (1: low, -1: high)  |
| ```time```    | Real          | 1         | time  |
| ```geomdata``` | GeometryData | -         |domain geometry data  |
| ```prob_parm```| ProbParm     | -         |problem parameter data  |


TO CLARIFY
```s_refl```


For example, a slip wall in y bottom wall woiudl be could be define as

```cpp
    if (idir == 1 && sgn == 1) {
        s_ext[URHO] = s_int[URHO];
        s_ext[UMX]  = s_int[UMX];
        s_ext[UMY]  = -s_int[UMY];
        s_ext[UMZ]  = s_int[UMZ];
        s_ext[UET]  = s_int[UET];
    }
```

If the boundary is not defined, the code will use the one specified in the input for that particular face.


# Immersed Boundaries

In **Immersed Boundary Methods** (IBM), boundary conditions on the surface are enforced by interpolating the influence of the boundary onto the fluid flow. This is achieved by assigning specific values at *ghost points* (GP). These values are extrapolated from the surrounding fluid points, taking into account the surface normal and the imposed boundary conditions.

<img src="../../images/IB_geom.png">


## Set-up 

### input

To set-up the Define STL file with geometry

```
#-------------------------------- IMMERSED BOUNDARY ---------------------------#
ib.filename = sphere_fine_20k.stl
```


### prob.h


```cpp
pp.add   ("ib.move",0); // 0=false, 1=true
pp.add   ("ib.plot_surf",0); // 0=false, 1=true
```

### tagging

### bcnormal



## Implementation: IBM class

```eib_t``` is explicit geometry (triangulation based) immersed boundary method
class. It holds an array of IBMultiFab, one for each AMR level; and it also holds
the geometry

### Markers

The field variable ibMarkers (is a boolean)  has two values 
```ibMarkers(i,j,k,0)``` is the **SLD** solid point marker and
```ibMarkers(i,j,k,1)```  the **GP** ghost  point marker. Both markers indicated if the i,j,k cell in teh domain is a solid or a ghost point.Remember that **GP** cells *are* **SLD** cells.
