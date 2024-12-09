# Thermo Solver

## Template Thermodynamics

Thermodynamics template is quite large, usually only requiring the index of the variables. It is usually built with some public parameters (gamma,mw, etc.) and a set of functions that are reuired by numerical solvers

```cpp
template <typename idx_t>
class calorifically_perfect_gas_t {
 protected:
 public:
  Real gamma = 1.40;   // ratio of specific heats
  Real mw = 28.96e-3;  // mean molecular weight air kg/mol
  Real Ru = gas_constant;
  Real cv = Ru / (mw * (gamma - Real(1.0)));
  Real cp = gamma * cv;
  Real Rspec = Ru / mw;

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void RYP2E(const Real R,
                                                      const Real* /*Y*/,
                                                      const Real P,
                                                      Real& E) const {
    E = P / (R * (gamma - Real(1.0)));
  }
```

### Functions

| **Function**     | **Feature**                                                             |
| ---------------- | ----------------------------------------------------------------------- |
| `RYP2E`          | computes specific internal energy as a function of pressure and density |
| `RYE2TP`         | computes temperature and pressure as a function of density and energy   |
| `RYE2Cs`         | computes sound speed                                                    |
| `cons2eigenvals` | compute eigenvalues array from array of conervative variable            |
| `prims2fluxes`   | compute fluxes from primtive array                                      |
| `cons2prims`     | converts conservative array to primitive variables array                |
| `prims2cons`     | converts primitve to conservative variables                             |
| `max_char_speed` | calcualtes maximum charatersic wave speed                               |
| `cons2char`      | converts conservative array to characteristic variables array           |
| `char2cons`      | converts characteristic variable array to consetvative variables array  |

The calls of this functions are (all declared `AMREX_GPU_DEVICE AMREX_FORCE_INLINE`)

```cpp
void RYP2E(const Real R, const Real* /*Y*/, const Real P, Real& E) const {
```

```cpp
void RYE2TP(const Real R,const Real* /*Y*/,const Real E, Real& T, Real& P) const {
```

```cpp
void RYE2Cs(const Real /*R*/,const Real* /*Y*/,const Real E,Real& cs) const {
```

```cpp
GpuArray<Real,idx_t::NWAVES> cons2eigenvals(const int i, const int j, const int k, const Array4<Real>& cons, const GpuArray<int, 3>& vdir) const {}
```

```cpp
void prims2fluxes( int i, int j, int k, const Array4<Real>& prims, const Array4<Real>& fluxes, const GpuArray<int, 3>& vdir) const {
```

```cpp
void cons2prims(const MFIter& mfi, const Array4<Real>& cons, const Array4<Real>& prims) const {
```

```cpp
void prims2cons(const IntVect& iv, const Array4<Real>& prims,Real cons[idx_t::NCONS]) const {
```

```cpp
Real max_char_speed(const IntVect& iv, const int dir, const int ng,
                    const Array4<const Real>& prims) const {
```

```cpp
void cons2char(RoeAvgState r, Real f[idx_t::NCONS]) const {
```

```cpp
void char2cons(RoeAvgState r, Real f[idx_t::NCONS]) const {
```

## Template Transport

Transport properties are built similarly, for example a constant conductivity

```cpp
template <typename param>
class cond_const_t {
 private:
 public:
  Real cond_ref = param::conductivity;
  
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real cond(Real& T) const {
    return cond_ref;
  }
};
```

This template requires an input parameter (as a **struct**) that defiens the parameter. Viscous solver will use a `cond(T)` function
