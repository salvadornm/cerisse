# LES models

This page explains how basic Large Eddy Simulation (LES) models are implemented in the `LES_t` class, defined in the`LES.h` file. The typical set-up is in `prob.h`

```cpp
typedef closures_dt<indicies_t, visc_suth_t, cond_suth_t,
                    calorifically_perfect_gas_t<indicies_t>, Smagorinsky_t<LESparm_t,indicies_t> > ProbClosures;
```

Where Smagorinsky is a LES-type template

```cpp
template <typename param, typename cls_t>
class LES_t {
```

The parameter structure defines constants and options for the model, for example

```cpp
  Real Prsgs = param::Prsgs;
  Real Scsgs = param::Scsgs;
  ...
```

These parameters include sub-grid Prandtl number , Schmidt number  (see [theory](../equations.md)).\
An example of a `LESparm_t` in `prob.h` is:

```cpp
struct LESparm_t {

  public:

  static constexpr int  order = 2;       // order numerical scheme (for gradient estimation)
  static constexpr Real Prsgs = 0.7;     // sgs Pr
  static constexpr Real Scsgs = 0.7;     // sgs Sc
  static constexpr Real Cs = 0.1;        // Smagorinsky constant
  static constexpr Real CI = 0.08;       // Yoshizawa constant
  static constexpr bool fixDelta = true; // Fix filter witdth
  static constexpr Real Delta = 0.02;    // Filter width (if above true)  L/20
};
```

The `defaultparm_t` option can still be selected by includinng `NumParam.h` (see default ioptions there)

The LES-type template must contain several functions:

```cpp
Real calc_delta( ...); // In LES 
void visc_sgs( ... , Real& mu_T ); //  In Smagorinsky/WALE/ ... model dependent
void cond_sgs( ... , Real& cond_T );
void diff_sgs( ... , Real& rhoD_T );
void compute_sgsterms( ... ,  Real& mu_T, Real& cond_T, Real& rhoD_T );
void void tau_sgs(... , Real& tau_T );
```

That provide the relevant sub-grid transport paramaters due to small scale mixing.\
These functions are called from `viscous.h`

## Smagorinsky

In the simplest Smagorinsky model, the sub-grid viscosity function is

```cpp
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void visc_sgs(
  const int i, const int j, const int k, const Array4<const Real>& q,
  const GpuArray<Real, AMREX_SPACEDIM>& dxinv, const Real delta, Real& mu_T)
```

This start by computing the velocity gradient tensor using central differences, given the LES order.

```cpp
    const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
    Real dUdx[3][3] = {{0.0}};
    // finite difference central order 2/4/6
    for (int m = 0; m < AMREX_SPACEDIM; m++) {
      for (int n = 0; n < AMREX_SPACEDIM; n++) {
        dUdx[m][n]  = normal_diff<param::order>(iv, n, this->QUn[m], q, dxinv); // dUmdn
     }
    }
```

The default number of ghost points in Cerisse is three , which limits order in standard central difference\
schemes to _sixth_ order. If the sub-grid viscosity is to be calculated in the ghost points, for example to be added to the viscosity the order of viscous term and LES cannot be independent.\
For example, if diffusive fluxes are to be solved with fourth order. The flux at the face is

$$
f_{i+1/2} = \left. \rho D \frac{\partial \phi}{\partial x} \right|_{i+1/2} \approx
\left( \alpha_{2} [\rho D]_{i+2} +  ... +  \alpha_{-2} [\rho D]_{i-2} \right) 
\frac{\beta_{2}\phi_{i+2} + ... + \beta_{-2}\phi_{i-2}}{\Delta x}
$$

where $$\alpha$$ and $$\beta$$ are interpolation and derivative face-centred weights respectively, see Section [numerics](../theory/equations/numerical-methods.md#central-differences).\
However, the sub-grid viscosity $$\left. \rho D \right|_{i+2}$$, requires evaluation of the velocity gradients

$$
\left. \rho D \right|_{i+2} = \mathcal{F} \left( \beta_{l} u_{i+2+l} + ... \beta{-l} u_{i+2+l} \right)
$$

Where the maximum $$l$$ depends on the LES order (1 for second order, 2 for fourth order, etc.).\
So $$i+2+l < i+NG$$, where $$NG$$ the number of ghost points, which by default is 3.\
In general

$$
\mbox{half stencil viscous}  + \mbox{half stencil LES}  \leq 3
$$

Therefore, _without changing the default value of GHOST POINTS(3)_

* For order 2 evaluation of viscous derivatives, LES order is limited to 2 or 4
* For order 4 evaluation of viscous derivatives, LES order is limited to 2
* For order 6 evaluation of viscous derivatives, LES cannot be used

{% hint style="warning" %}
The default number of GHOST POINTS can be by changing `indicies_t` in the `prob.h` to the more general`indices_gent<>` passing the number of GHOST points as argument. Beware using this, as the number of GHOST points is used in many places, including I/O subroutines.
{% endhint %}
