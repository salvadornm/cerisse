# Code Structure



These pages are for users who want to tinker with the numerical solvers and the code structures to create
new thermodynamics, numerical methiods.

Although a deep understanding of the AMReX structure is not necessary, it is recommended to familiarize yourself with how AMReX handles data in general. 

The code extensively utilizes templates, sometimes in unconventional ways. 

Initially, there is a basic template that defines, how many scalars are to be solved and
where are located in the storage array. It defines a lot of indicators of where each variable is. 
This is done in **indicies_t**.

This template is combiend with templates for the thermodynamics, transport properties and turbulence models (if present) define the **closures_dt** template. This is done in `prob.h`, where we define a
**ProbClosures** 

```cpp
typedef closures_dt<indicies_t, visc_suth_t, cond_suth_t,
                    calorifically_perfect_gas_t<indicies_t>> ProbClosures;
```
or

```cpp
using ProbClosures = closures_dt<indicies_stat_t, visc_suth_t, cond_suth_t,
                                 calorifically_perfect_gas_t<indicies_t>>;
```

Each of this templates corresponds to a file (most of the times)
For example the Sutherland viscosity model defined in `src/cls/Transport.h`

```cpp
class visc_suth_t {
 private:
  // Sutherland's fit from Computational Fluid Mechanics and Heat Transfer
  Real visc_ref = 1.458e-6;
  Real Tvisc_ref = 110.4;

 public:
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real visc(const Real& T) const {
    return visc_ref * T * sqrt(T) / (Tvisc_ref + T);
  }
};
```