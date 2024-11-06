# PROB


The file is structured in several parts, not all of them are necessary.
Use ```prob.h``` form one of the examples.




## Problem parameters

```cpp
struct ProbParm {
  Real p_l = 1.0;
  Real p_r = 0.1;
  Real rho_l = 1.0;
  Real rho_r = 0.125;
  Real u_l = 0.0;
  Real u_r = 0.0;
};
```

### Global Cosntants

Similarly,  global constants can be defined

```cpp
static constexpr Real Reynolds = 3000.0;  
```

Define them as ```static constexpr```, so 
they are evaluated at compile time, making it efficient for use in other compile-time expressions.
For example, they may be used to set-up a constant viscosity


## Themodynamic and Transport Closures

creating the "closures" of the problems, that means which model of transport 
properties is to be used. 

For example:

```cpp
typedef closures_dt<indicies_t, visc_suth_t, cond_suth_t,
                    calorifically_perfect_gas_t<indicies_t>> ProbClosures;
```

This selects a problem that uses the structure `indicies_t`, which tell position
of the variables and how many variables to solve. The above also selects
Sutherland viscosity model  ```visc_suth_t```
All this is wraped in the class ```ProbClosures``


Options in  ***closures_dt***


| Closure                      | Type          | Use | Description                                                  |
| --------------------------- | ------------- |:-------:| ------------------------------------------------------------ |
| ```visc_suth```             |            |   often      | Sutherland Viscosity                           |
| ```cond_suth```             |            |   often      | Sutherland Conductivty                         |
| ```calorifically_perfect_gas_t```             |            |   often      | perfect gas                      |


### Passing Arguments

Parameters can be passed by a small structure `methodparm_t`,
so all paraeters  can be changed directly in `prob.h`, 
such as the order of the scheme, disspation parameters, 
extra damping, sensor variables, etc. The above example would be

```cpp
struct methodparm_t {

  public:

  static constexpr int  order = 2;              // order numerical scheme viscous
  static constexpr Real conductivity = 0.0262;  // conductivity (for constant value)
  static constexpr Real viscosity   = 1.85e-5;  // viscosity    (for constant value)
  
};
typedef closures_dt<indicies_t, visc_const_t<methodparm_t>, cond_const_t<methodparm_t>,
                    calorifically_perfect_gas_t<indicies_t>> ProbClosures;
```


## Equations

Define the RHS of the problem, that includes which equation to solve and which numerivcal scheme to use. In general 

$$
\frac{\partial U}{\partial t} = \mbox{RHS}(U)
$$

where the RHS includes, inviscid (Euler) terms viscous (for Navier-Stokes) and source terms. For example:

```cpp
typedef rhs_dt<skew_t<true,false, 4, ProbClosures>, no_diffusive_t, no_source_t > ProbRHS;
```
This will solve
$$
\frac{\partial U}{\partial t} = F(U)
$$
corresponding to the Euelr equations, with a 4th order skew-symmetric numerical scheme.
Template `skew_t`,`diffusive_t` may have  arguments (like the order of the scheme).




This call indicates the different type of numerical scheme for the advection part, default

```cpp
  typedef rhs_dt<keep_euler_t<false, false, 4, ProbClosures>,no_diffusive_t,
               no_source_t> ProbRHS;
```

| RHS                     | Options          | Use | Description                                                  |
| --------------------------- | ------------- |:-------:| ------------------------------------------------------------ |
| ```riemann_t```             |            |   often      | Riemann Solver   MUSCL              |
| ```keep_euler_t```             |  AD,IB, order           |   often      | KEEP Scheme                 |
| ```centraldif_t```             |  AD,IB, order           |   often      | Central Scheme                 |
| ```skew_t```             |  AD, IB, order           |   often      | Skew-symmetric Scheme                 |
| ```no_diffusive_t```        |            |   often      | No diffusive part (Euler)                |
| ```no_source_t```           |            |   often      | No source term    |



## Initial Condition

The **prob_initdata** function, this is function is called for every 
```i,j,k``` cell

```cpp
// initial condition
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_initdata(int i, int j, int k, Array4<Real> const &state,
              GeometryData const &geomdata, ProbClosures const &cls,
              ProbParm const &prob_parm) {
```

In most examples,  a few auxiliar definitions follow that extract the size the domain
and allow to define the spatial coordinate *x* (and *y* and *z*) of the cell

```cpp
  const Real *prob_lo = geomdata.ProbLo();
  const Real *prob_hi = geomdata.ProbHi();
  const Real *dx = geomdata.CellSize();
  Real x = prob_lo[0] + (i + Real(0.5)) * dx[0];
```

NOTE: Real for variable precision (check ```GNU_Makefile```)


This allow to define, for example, different condtion depending on the  *x* coordinate.
The structure ```prob_parm```of type  **ProbParm** is used to recall problem parameters


```cpp
  Real Pt, rhot, uxt;
  if (x < prob_hi[0] / 2) {
    Pt = prob_parm.p_l;
    rhot = prob_parm.rho_l;
    uxt = prob_parm.u_l;
  } else {
    Pt = prob_parm.p_r;
    rhot = prob_parm.rho_r;
    uxt = prob_parm.u_r;
  }
  Real et = Pt / (cls.gamma - Real(1.0));
```

The function needs to fille the ```state``` array, where the conservative variables exist.
The index **URHO**, **UMX** are defined in the template ```indicies_t``` and are accessed through the class
as ```cls.URHO``` and so on.

``` cpp
  state(i, j, k, cls.URHO) = rhot;
  state(i, j, k, cls.UMX)  = rhot * uxt;
  state(i, j, k, cls.UMY)  = Real(0.0);
  state(i, j, k, cls.UMZ)  = Real(0.0);
  state(i, j, k, cls.UET)  = et + Real(0.5) * rhot * uxt * uxt;
}
```

Since this is pure C++ code, it allows for the construction of complex initial conditions. PelePhysics can also be utilized to incorporate chemistry profiles, such as those for 1D premixed flames.

NOTE: indicies is common, an dit allocates spaces for so in 1D  ```state(i,j,k,cls.UMY)=0```, 
```state(i,j,k,cls.UMZ)=0``` and in 2D UMZ=0




## Source

Source term


## Boundary Condition

The ```prob.h``` can be used to implement user-specific boundary conditions (transient, turbulent, etc).
To do that,  a local template is defined



## Other

The names follow in the **cons_var_names** array

```cpp
inline Vector<std::string> cons_vars_names={"Xmom","Ymom","Zmom","Energy","Energy"};
```

The type of variables, keep as it is, scalar ser 0 and vectors are given by their components


```cpp
inline Vector<int> cons_vars_type={1,2,3,0,0};
```


Data missed from input file

```cpp
void inline inputs() {
  ParmParse pp;

  pp.add("cns.order_rk", 3);   // -2, 1, 2 or 3"
  pp.add("cns.stages_rk", 3);  // 1, 2 or 3
}
```
