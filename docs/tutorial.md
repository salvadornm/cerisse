# Tutorial

This  page will explain the set-up, run and visualization of a typical case.
For a quick overview of running the code, see [Quickrun](quickstart.md).
At a glance, domain dimension and control parameters are handled in file `inputs` (or similar name)  while problem description is in * `prob.h`

NOTE: It is recommended to cp the tutorial into a `wrk` or `exm/tmp` directory. Is not neccesary, but if the code is cloned, changes in mian files will register as to commit.
The previous directories will always be ignored by **git**.


## Problem Set-up

The problem involves a heavy fluid descending into a lighter fluid, creating conditions for a Rayleigh-Taylor instability. The upper half of the domain contains a fluid with a density of 2, while the lower half is filled with a fluid of density 1. The initial pressure distribution is hydrostatic, and a velocity perturbation is introduced to trigger the instability
The initial conditions can be summarised as

$$
 (\rho, v_x, v_y, P) =\left \{
  \begin{array}[cr]  
  ( 2, 0, -\epsilon c \cos(8 \pi x),P_0 + \rho g y )   && \mbox{ at }  y > y_0  \\
  ( 1, 0, -\epsilon c \cos(8 \pi x),P_0 + \rho g y )   && \mbox{ otherwise }   \\
  \end{array}
\right.
$$

where Inline equation: *y0* is the middle of the domain, *P0=2* gravity is taken as 1 pointing downwards and the perturbation is 0.025.
The initial conditions followed the paper by Shi et al [^1].
Both fluids are miscible and follow the ideal gas law with adiabatic coefficient of 5/3.


## prob.h


### Set-up of Problem Parameters

The parameters of the problem (not the domain) are wrapped into a **ProbParm** structure

```cpp
// problem parameters
struct ProbParm {
  Real p_int = 2.0;
  Real rho_1 = 1.0;
  Real rho_2 = 2.0;
  Real grav = -1.0; 
  Real eps =  0.025;
};
```

### Set-up of  initial conditions

This is done in the function **prob_initdata**

```cpp
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_initdata(int i, int j, int k, Array4<Real> const &state,
              GeometryData const &geomdata, ProbClosures const &cls,
              ProbParm const &prob_parm) {
  const Real *prob_lo = geomdata.ProbLo();
  const Real *prob_hi = geomdata.ProbHi();
  const Real *dx = geomdata.CellSize();

  Real x = prob_lo[0] + (i + Real(0.5)) * dx[0];
  Real y = prob_lo[1] + (j + Real(0.5)) * dx[1];
```
where the spatial coordinates **x** and **y** of the cells are determined from the mesh sizes and problem dimensions (which are defined in **inputs**).

The initial conditions are defined as

```cpp
  Real Pt, rhot, uxt,uyt;
  Real Lint = prob_hi[0] / 2;
  Real Pint = prob_parm.p_int; // interface Pressure

  const Real freq = Real(8)*Real(3.14159265359); // wavelength = x-domain

  Real yrel = y - Lint;
  Real delta= 0.2*Lint;   // region size where perturbation is significant
  Real delta2  = dx[1]/5; // transition region between top/bottom
  Real step = Real(0.5) + Real(0.5)*tanh(yrel/delta2);
  rhot = step*prob_parm.rho_2 + (Real(1.0) -step)*prob_parm.rho_1;
  Pt = Pint + rhot*prob_parm.grav*(y - Lint); // hydrostatic pressure

  uxt = Real(0.0);
  uyt = -prob_parm.eps*cos(freq*x)*aux; // perturbation in y-component

  state(i, j, k, cls.URHO) = rhot;
  state(i, j, k, cls.UMX) = rhot * uxt;
  state(i, j, k, cls.UMY) = rhot * uyt;
  Real et = Pt / (cls.gamma - Real(1.0));
  state(i, j, k, cls.UET) = et + Real(0.5) * rhot * (uxt * uxt + uyt * uyt); 
```


The array **state** holds all the information required in conserved variables

$$
\mbox{state} = \left( \rho, \rho v_x, \rho v_y, \rho e_T  \right)
$$

to initialise the problem.

### Solvers

The solvers is selected in the line

```cpp
typedef rhs_dt<riemann_t<false, ProbClosures>, no_diffusive_t, user_source_t<ProbClosures>>
    ProbRHS;
```
which defines which problem to solve


### Source Term

Function

```cpp
template <typename cls_t>
class user_source_t {
  public:
  void inline src(const amrex::MFIter &mfi,
                  const amrex::Array4<const amrex::Real> &prims,
                  const amrex::Array4<amrex::Real> &rhs, const cls_t *cls_d,
                  amrex::Real dt){
```

The lines below loop over the cells and gravity is added to the momentum 
and energy equation. The `prob_parm.grav` reads the gravity value defined
in `ProbParm` structure. 

```cpp
amrex::ParallelFor(bx,
  [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {
    const auto& cls = *cls_d;
    Real rho = prims(i, j, k, cls.QRHO);
    rhs(i,j,k,cls.UMY) += prob_parm.grav*rho;
    rhs(i,j,k,cls.UET) += prob_parm.grav*rho*prims(i, j, k, cls.QV);
    });
```


### AMR

This function defines which points need to be "tagged" to be refined.
AMReX will then create patches around that cell to create the refined regions.

```cpp
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
user_tagging(int i, int j, int k, int nt_level, auto &tagfab,
             const auto &sdatafab, const auto &geomdata,
             const ProbParm &prob_parm, int level) {
```
The actual tagging is produced in the follwoing lines, where the
array `tagfab(i,j,k)` is set ot true or false based on the 
value of density.

```cpp
  Real rhot = sdatafab(i,j,k,ProbClosures::URHO);
  if (nt_level > 0)
  {
  // tag cells based on density values
  switch (level)
   {
    case 0:
      tagfab(i,j,k) = (rhot > 1.1 && rhot < 1.9);
      break;
    case 1:
      tagfab(i,j,k) = (rhot > 1.2 && rhot < 1.8);
      break;
    default:
      tagfab(i,j,k) = (rhot > 1.3 && rhot < 1.7);
      break;
    }
  }  
```

## Input file


### Time steps
```
max_step  = 20000
stop_time = 2.0
cfl = 0.3
```

### Boundary conditions

The boudantry conditions are defined
More complex boudnarty conditions can be defined 

```
cns.lo_bc = -1 5
cns.hi_bc = -1 5
```

### Mesh and refinement

```
amr.n_cell          = 128 512
amr.max_level       = 1     # maximum level number allowed
amr.regrid_int      = 2 2 # how often to regrid per level
```
These lines will select size of the mesh and how often to regrid



## Compile and Running

To compile use

```bash
$ make
```

and run case using 2 cores (or any number you want)

```bash
$ mpirun -np 2 ./main2d.gnu.MPI.ex inputs
```

## Post-processing

The input file is set up to plot every 100 steps

```bash
$ ls plot
plt00000	plt01300	plt02600	plt03900	plt05200
plt00100	plt01400	plt02700	plt04000	plt05300
plt00200	plt01500	plt02800	plt04100	plt05400
plt00300	plt01600	plt02900	plt04200	plt05500
plt00400	plt01700	plt03000	plt04300	plt05600
plt00500	plt01800	plt03100	plt04400	plt05700
plt00600	plt01900	plt03200	plt04500	plt05800
plt00700	plt02000	plt03300	plt04600	plt05900
plt00800	plt02100	plt03400	plt04700	plt06000
plt00900	plt02200	plt03500	plt04800	plt06100
plt01000	plt02300	plt03600	plt04900	plt06200
plt01100	plt02400	plt03700	plt05000	plt06300
plt01200	plt02500	plt03800	plt05100	plt06390
```

In the following sectionsm we will show how to visualize the results using Python,
Visit and Paraview.

### Python

The following assumes that **yt** is installed. Check [Tips](tips.md) for installation.
A Python script for easy of use,by invoking 

```bash
$ python plot.py
```

The last density snaphot at *t=2* will be saved into a png file, which looks like
(with default color scheme in **yt**).

<img src="../../images/tutorial_yt1.png" width=50%>


### Visit

To see the results with [VisIt](https://visit-dav.github.io/visit-website/),  you can use the script **cerisse**, to open all directories at the same time (to make an animation for example). See  [Tips](tips.md) to set-up the script.

```bash
$ cerisse visit
```
This will create a file `movie.visit` in the tutorial directory.
To open Visit, type `visit` (or use the appropiate icon).

Once opened, use **File/Open** to load the  `movie.visit file.
If it loads corerctly, it should look something like:

<img src="../../images/tutorial_Visit1.png" width=100%>


Adding a plot is easy by just  **Add/Pseudocolor/Density**

<img src="../../images/tutorial_Visit2.png" width=100%>

The results shows the pseudocolor for density with default visualization options
(see Visit manual and webpage for details on costumization).
By clicking the arrow in the panel
<img src="../../images/tutorial_Visit3.png" width=50%>

we obtain  a quick animation of the results. The final snapshot
should be

<img src="../../images/tutorial_Visit4.png" width=100%>

By **Add/Subset/Levels** we can get an idea of the refinement 

<img src="../../images/tutorial_Visit5.png" width=50%>

Double-cliking on the **Subset/Levels**  we can edit the attributes of the plot

<img src="../../images/tutorial_Visit6.png" width=50%>

And the final plot  should look like, where the solid lines are the limits of refinmement

<img src="../../images/tutorial_Visit7.png" width=100%>

There is no need to open the complete database and tiem steps can be opened individually by open the individual header files (for example opening directly `plot/plt02600/Header`)
VisIt is good for exploring the data interactively, Visit can also be scripted with Python.
For a VisIt tutorial check (http://visitusers.org/index.php?title=VisIt_Tutorial)


### Paraview

To use [Paraview](https://www.paraview.org), start Paraview in the usual way and open the `plot\plt...` folder


<img src="../../images/tutorial_Paraview1.png" width=100%>

Select the AMReX/Boxlib Grid Reader 

<img src="../../images/tutorial_Paraview2.png" width=50%>

Tick the density box and click **Apply**, it would look something like

<img src="../../images/tutorial_Paraview4.png" width=100%>

Select Density and Surface

<img src="../../images/tutorial_Paraview3.png" width=100%>


And then use the  final time to show the final plot (with default options)

<img src="../../images/tutorial_Paraview5.png" width=100%>

Similar to Visit, Paraview is good for exploring the data interactively.
can also be scripted with Python. Both softwares are similar and can 
use HPC, remote visualization and support large number of points (billions)
A Paraview Tutorial manual can be found in [Paraview Manual](https://www.paraview.org/Wiki/images/5/5d/ParaViewTutorial41.pdf)




[^1]:  J Shi et al, *Resolution of high order WENO schemes for complicated flow structures*, J Computational Physics, (2003), 186, pp 690-696.
https://www.sciencedirect.com/science/article/pii/S0021999103000949
