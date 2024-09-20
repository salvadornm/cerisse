# Classes


# rhs/

In the prob.h file, there the **rhs_dt** is defined

```cpp
typedef rhs_dt<riemann_t<false, ProbClosures>, no_diffusive_t, no_source_t>
    ProbRHS;
```
The line above  defines the type of problem to be solved, and consist of a serie of keywords that represent the relative class


| class                      | Type          | Options| Description                                                  |
| --------------------------- | ------------- |:-------:| ------------------------------------------------------------ |
| riemann_t              | Euler Solver           |   bool, ProbClosures     | Class for Riemann Solvers   |
| rusanov_t              | Euler Solver           |   bool, ProbClosures     | Rusanov   |
| keep_t              | Euler Solver           |   bool, ProbClosures     | KEEP scheme   |
| skew_t              | Euler Solver           |   bool, ProbClosures     | Skew-symmetric scheme   |
| no_diffusive_t             | Option          | 0       | No diffusion, viscosity   |
| no_source_t             | Option          | 0       | No source term   |
| user_source_t             | Option          | 0       | Defined in prob.h    |




## Riemann

Template class for Riemann Solver, in file **Riemann.h**

```void riemann_prob```

Solves the HLLC  Riemann Problem based on left and right states


### Limiter 

```Real limiter(dlft,drgt)```

MUSCL type reconstruction with Vanleer limiter

```cpp
  Real dcen = Real(0.5) * (dlft + drgt);
  Real dsgn = Math::copysign(Real(1.0), dcen);
  Real slop = Real(2.0) * min(Math::abs(dlft), Math::abs(drgt));
  Real dlim = (dlft * drgt >= Real(0.0)) ? slop : Real(0.0);
  return dsgn * min(dlim, Math::abs(dcen));
```

The average gradient is defined by

$$
d_{cen}= \frac{1}{2} ( d_{left} + d_{right})
$$

Using the minmod limiter, only applies if $$ d_{left}d_{right} >0 $$
(no extrema, otherwise 0)

$$
d_{lim} = 2  \min(|d_{left}|,|d_{right}|)
$$


### Reconstruction

```void cns_slope_x(i,j,k,&dq,&q,&cls)```


Example of reconstruction in X, using characteristics.
This fucntion calculates the limited gradients 

For example for u-x xcompoenent


```cpp
dlft =  Real(0.5) * (q(i, j, k, cls.QPRES) - q(i - 1, j, k, cls.QPRES)) /
               cspeed +
        Real(0.5) * q(i, j, k, cls.QRHO) *
               (q(i, j, k, cls.QU) - q(i - 1, j, k, cls.QU));
drgt =  Real(0.5) * (q(i + 1, j, k, cls.QPRES) - q(i, j, k, cls.QPRES)) /
               cspeed +
        Real(0.5) * q(i, j, k, cls.QRHO) *
               (q(i + 1, j, k, cls.QU) - q(i, j, k, cls.QU));
Real d2 = limiter(dlft, drgt);
```
which corresponds to 

$$
d_{left} = \frac{1}{2}  \frac{ p_i - p_{i-1} }{c}  +   \frac{1}{2} \rho_i \left(u_i - u_{i-1} \right) 
$$
$$
d_{right} = \frac{1}{2} \frac{ p_{i+1} - p_i }{c}  +   \frac{1}{2} \rho_i \left(u_{i+1} - u_i \right) 
$$

The slope is then recointructed using the limiter

$$
d_2 = \Phi(d_{left},d_{right})
$$

same applies for y and and z.


### Solving Riemann problem x

```void cns_riemann_x```

-Reconstruct using characteristics

-Return to primitives

-Call Riemann solver


## Central


## Weno


## CentralKEEP




