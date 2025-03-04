# Viscous and Diffusion solver

This page explains how viscosity and diffusion are implemented in Cerisse 
in the ```viscous_t``` class, defined in  the```viscous.h``` file

The main objective of this class is to provide
the function `dflux()` (or `dflux_ibm()` for immersed boundary methods) to
calculates **diffusive fluxes** that are added to the flux array.

```cpp
  void inline dflux(const Geometry& geom, const MFIter& mfi,
            const Array4<Real>& prims, std::array<FArrayBox*, AMREX_SPACEDIM> const &flxt, 
            const Array4<Real>& rhs, const cls_t* cls) {
```

## Summary

### **1. Class Template and Initialization**
- The class is templated on `param` and `cls_t`, allowing flexibility for different numerical orders and physics closures.
- It initializes numerical coefficients (`CDcoef`, `INTcoef`) for discretization schemes.
- The constructor `viscous_t()` calls `calc_CDcoeffs<param::order>()` to compute these coefficients.

### **2. Transport Property Calculation**
- It determines transport properties such as:
  - **Dynamic viscosity** (`mu_arr`)
  - **Thermal conductivity** (`lam_arr`)
  - **Bulk viscosity** (`xi_arr`)
  - **Species diffusivity** (`rhoD_arr`)
- If `PelePhysics` is enabled, transport coefficients are obtained using **PelePhysics Transport Models**; otherwise, they are computed using `cls->visc()` and `cls->cond()`.

### **3. Diffusion Flux Computation**
- `cns_diff()` computes **viscous, heat, and mass diffusion fluxes** using:
  - Finite difference approximations (`normal_diff`, `tangent_diff`).
  - Interpolation functions (`interp<param::order>()`).
  - Divergence of velocity (`divu`) and viscous stress components (`tau11`, `tau12`, etc.).
  - Energy flux contributions from species diffusion.
- If immersed boundaries (`IBM`) or embedded boundaries (`EB`) are present, `cns_diff_ibm()` modifies the approach by reducing accuracy near boundaries to maintain numerical stability.

## Diffusion flux implementation

The diffusion flux follows

$$
\vec{J}_k = \rho Y_k \vec{V}_k 
$$

where

$$
\vec{V}_k = D_k \left[ \vec{\nabla} X_k +  \left(X_k - Y_k \right) \vec{\nabla}  \ln p \right]
$$

The implementation of Transport properties in [PelePhsyics](https://pelephysics.readthedocs.io/en/latest/Transport.html) is based on 
[ERG](https://doi.org/10.1006/jcph.1995.1151) to compute mixture-averaged diffusivities for each species. All this is handled by PelePhysics through FUEGO model defined at compile time.
There is no species tarnsport if only one species present.


### Correction to maintain mass

By mass conservation, the sum of mass fluxes must be 0 
$$
\sum_k  \vec{J}_k = 0 
$$

However, there is no guarantee of this. Following [PeleC](https://amrex-combustion.github.io/PeleC/Algorithms.html#diffusion)  arbitrary correction flux is added to all species

$$
\vec{V}_c = \sum_k Y_k \vec{V}_k 
$$
The new *corrected* flux operates with a modified velocity
$$
\vec{J}^{corr}_k = \vec{J}_k -\rho Y_k \vec{V}_c  = \rho Y_k \vec{V}_k  - \rho Y_k V_c
$$
and
$$
\sum_k \vec{J}^{corr}_k = \rho \underbrace{ \sum_k  Y_k \vec{V}_k}_{V_c} - \rho V_c \underbrace{\sum_k Y_k}_{1} = 0
$$

The energy flux is also affected by energy transport due to mass transport (due to chemical species enthalpy)

The actual code lines are (mole, mass and enthakopy arrays computed before)

```cpp
    const Real dpdx  = normal_diff<param::order>(iv, d1, cls_t::QPRES, q, dxinv); 
    const Real pface = interp<param::order>(iv, d1, cls_t::QPRES, q);
    const Real dlnp = dpdx/pface;
     
    Real Vc = 0.0;
    for (int n = 0; n < NUM_SPECIES; n++) {
      Real Xface = 0.0,Yface= 0.0,hface=0.0,dXdx=0.0;
      for (int l = 0; l < param::order; l++) {
        Xface += xmole[l][n]*INTcoef(l);
        Yface += ymass[l][n]*INTcoef(l);
        hface += hi[l][n]*INTcoef(l);
        dXdx  += xmole[l][n]*CDcoef(l);
      }      
      dXdx  /= dxinv[d1];
      const Real rhoD_f = interp<param::order>(iv, d1, cls_t::CRHOD + n, coeffs); 

      const Real Vd = -rhoD_f * (dXdx + (Xface - Yface) * dlnp);
      Vc += Vd;
      flx(iv, cls_t::UFS + n) += Vd; 
      flx(iv, cls_t::UET)     += Vd * hface;
    }
    // Add correction velocity to fluxes so sum(Vd) = 0
    for (int n = 0; n < NUM_SPECIES; ++n) {
      flx(iv, cls_t::UFS + n)-= Yface * Vc;
      flx(iv, cls_t::UET)    -= Yface * hface * Vc;
    }

```