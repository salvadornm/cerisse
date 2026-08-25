# NSCBC species convection

This two-dimensional PelePhysics case follows the species-convection test in
the supplied NSCBC+AMReX paper. The domain contains a Gaussian oxygen mass
fraction,

```text
Y_O2 = exp(-[(x-a)^2+(y-b)^2]/[2 R^2]),  Y_N2 = 1-Y_O2,
```

at constant pressure and temperature. It is convected toward the x-high NSCBC
outlet at 20 m/s. The x-low type-1 NSCBC supplies dry air
(`Y_O2=0.233`, `Y_N2=0.767`), while the type-3 x-high boundary lets the oxygen
distribution leave with pressure relaxation. The y direction is periodic to
isolate species transport from wall boundary layers.

The case uses the two-species PelePhysics `air` mechanism and the fourth-order
skew-symmetric convective scheme. It has no reactions or molecular diffusion.
Alternatively a KEEP-style 4th order can be used (uncomment the line in prob.h) 

