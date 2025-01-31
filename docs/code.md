# Overview

These pages are intended for users interested in exploring and customizing the numerical solvers and code structures to develop new thermodynamic models or numerical methods. While a deep understanding of the AMReX framework is not strictly required, it is highly recommended to familiarize yourself with how AMReX manages data structures and operations.

The code makes extensive use of templates, often employing them in unconventional ways. At its core, a basic template is used to define the number of scalar variables to be solved and their respective locations in the storage array. This template specifies numerous indicators for the position of each variable and is encapsulated in the **indicies\_t** structure.

This foundational template is then combined with additional templates for thermodynamics, transport properties, and turbulence models (if applicable), which collectively form the **closures\_dt** template.

This is done in `prob.h`, where we define a **ProbClosures**

```cpp
typedef closures_dt<indicies_t, visc_suth_t, cond_suth_t,
                    calorifically_perfect_gas_t<indicies_t>> ProbClosures;
```

or

```cpp
using ProbClosures = closures_dt<indicies_stat_t, visc_suth_t, cond_suth_t,
                                 calorifically_perfect_gas_t<indicies_t>>;
```

Each of this templates corresponds to a file (most of the times) For example the Sutherland viscosity model defined in `src/cls/Transport.h`

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

## Code Development

Use AmReX methodology

To add a new feature, the procedure is:

1. Create a branch for the new feature (locally):

```bash
$ git checkout -b AmazingNewFeature
```

2. Develop the feature, merging changes often from the **develop** branch into your AmazingNewFeature branch:

```bash
$ git commit -m "Developed AmazingNewFeature"
$ git checkout develop
$ git pull                      # fix any identified conflicts between local and remote branches of "develop"
$ git checkout AmazingNewFeature
$ git rebase develop        # fix any identified conflicts between "develop" and "AmazingNewFeature" 
```

3. To mare sure the feature does not bertak the code do test

```bash
$ cerisse test
```
