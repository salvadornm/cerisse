#include "wall_model.H"

#include <AMReX_LUSolver.H>
#include <PelePhysics.H>

#ifndef WALLMODEL_DEBUG
#include "CNS.H"
#else
#include <PelePhysics.H>
#endif

using namespace amrex;

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void LawOfTheWall::parallel_wall_stress(
  Real u, Real T, Real rho, Real Y[NUM_SPECIES], Real h, Real mu, Real lam,
  Real T_wall, Real& tau, Real& q)
{
  // u^+ = 1/kappa * log(y^+) + B
  const Real kappa = 0.4;
  const Real B = 5.5;
  const int maxiter = 10;

  // assume mu = mu_w, rho = rho_w
  Real nu_w = mu / rho;

  // Newton-Raphson: f(ut) = u/ut - log(h*ut/nuw)/kappa - B = 0
  // ut^{n+1} = ut^n - f(ut^n)/f'(ut^n)
  Real ut = std::sqrt(nu_w * u / h); // initial guess u_tau
  Real dut = 1.e10;                  // residual of u_tau
  int iter = 0;
  while ((iter < maxiter) && (std::abs(dut) > u * 1.e-5)) {
    dut = (u / ut - std::log(h * ut / nu_w) / kappa - B) / // f(ut)
          (-u / ut / ut - 1.0 / kappa / ut);               // f'(ut)
    ut -= dut;
    iter++;
    // std::cout << iter << " " << ut << std::endl;
  }

  tau = rho * ut * ut; // wall shear stress
  if (T_wall < 0) {
    q = 0.0; // adiabatic wall
  } else {
    // T = Tw + (Tr - Tw)*u/ue + (Te - Tr)*(u/ue)^2 (Walz 1969)
    Real r = 0.9; // recovery factor (~Pr^1/3 ~0.9 for turbulent flow)
    Real cp;
    auto eos = pele::physics::PhysicsType::eos();
    eos.RTY2Cp(rho, T, Y, cp);
    Real Tr = T + 0.5 * r * u * u / cp;
    q = lam * (Tr - T_wall) / u * tau / mu;
  }
}

#ifdef WALLMODEL_DEBUG
namespace {
void printArray(Real* arr, int n, std::string name)
{
  std::cout << name << ": ";
  for (int i = 0; i < n; ++i) { std::cout << arr[i] << " "; }
  std::cout << std::endl;
}

template <int N>
void printArray(Array1D<Real, 0, N - 1> arr, std::string name)
{
  std::cout << name << ": ";
  for (int i = 0; i < N; ++i) { std::cout << arr(i) << " "; }
  std::cout << std::endl;
}

template <int N>
void printMatrix(Array2D<Real, 0, N - 1, 0, N - 1, Order::C> mat, std::string name)
{
  std::cout << name << ":\n";
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) { std::cout << mat(i, j) << " "; }
    std::cout << std::endl;
  }
}
} // namespace
#endif

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void EquilibriumODE::parallel_wall_stress(
  Real uwm, Real Twm, Real rhowm, Real Y[NUM_SPECIES], Real h, Real /*mu*/,
  Real /*lam*/, Real T_wall, Real& tau, Real& q)
{
  constexpr int n_grid = 30;
  const Real stretch = 1.1;
  const Real uw = 0.0; // wall velocity
  auto eos = pele::physics::PhysicsType::eos();
  auto trans = pele::physics::PhysicsType::transport();

#ifdef WALLMODEL_DEBUG
  pele::physics::transport::TransportParams<
    pele::physics::PhysicsType::transport_type>
    trans_parms;
  trans_parms.allocate();
  trans_parms.host_trans_parm().const_viscosity = 1.0e-6 * 1.225e-3;
  trans_parms.host_trans_parm().const_conductivity = 0.025;
  trans_parms.sync_to_device();
  auto const* ltransparm = trans_parms.device_trans_parm();
#else
  auto const* ltransparm = CNS::trans_parms.device_trans_parm();
#endif

  const Real kappa = 0.41; // von Karman constant
  const Real Aplus = 17.0; // constant in the log-law
  const Real Prt = 0.9;    // turbulent Prandtl number

  Real y[n_grid], u[n_grid], T[n_grid]; // solution arrays
  for (int i = 0; i < n_grid; ++i) {
    Real ylo = h * (std::pow(stretch, i) - 1) / (std::pow(stretch, n_grid) - 1);
    Real yhi = h * (std::pow(stretch, i + 1) - 1) / (std::pow(stretch, n_grid) - 1);
    y[i] = 0.5 * (ylo + yhi);
    u[i] = y[i] / h * (uwm - uw) + uw; // linear initial guess
    T[i] = T_wall > 0.0 ? y[i] / h * (Twm - T_wall) + T_wall : Twm;
  }
  // printArray(y, n_grid, "y");
  // printArray(u, n_grid, "u");
  // printArray(T, n_grid, "T");

  Real res_u = 1e10, res_T = 1e10; // residuals
  Real rho_0, tauw, tmp[n_grid], mu[n_grid], lam[n_grid], cp[n_grid], mut[n_grid]; // temporary variables
  int iter = 0, max_iter = 10;
  LUSolver<n_grid, Real> lusolver;
  while ((res_u > 0.01 * uwm || res_T > 0.01 * Twm) &&
         (iter < max_iter)) { // until within 1% error
    // compute coefficients
    for (int i = 0; i < n_grid; ++i) {
      Real rho_i = rhowm * Twm / T[i]; // assume const p
      Real get_mu, get_lam, dummy;
      trans.transport(false, true, true, false, false, T[i], rho_i, Y, nullptr,
                      nullptr, get_mu, dummy, get_lam, ltransparm);
      mu[i] = get_mu;
      lam[i] = get_lam;
      eos.RTY2Cp(rho_i, T[i], Y, cp[i]);
      if (i == 0) { 
        tauw = mu[0] * (u[0] - uw) / y[0]; 
        rho_0 = rho_i;
      }
      Real yplus_i = y[i] * std::sqrt(tauw / rho_0) / (mu[0] * rho_0);
      mut[i] = kappa * y[i] * std::sqrt(rho_i * tauw) *
               std::pow(1.0 - std::exp(-yplus_i / Aplus), 2);
    }
    // printArray(mu, n_grid, "mu");
    // printArray(yplus, n_grid, "y+");
    // printArray(mut, n_grid, "mut");
    // printArray(lam, n_grid, "lam");

    // solve for u
    Array2D<Real, 0, n_grid - 1, 0, n_grid - 1, Order::C> A;
    Array1D<Real, 0, n_grid - 1> b;
    for (int i = 0; i < n_grid; ++i) {
      for (int j = 0; j < n_grid; ++j) { A(i, j) = 0.0; }
      b(i) = 0.0;
    }
    // wall BC
    Real mu_lo = mu[0] + mut[0];
    Real mu_hi = 0.5 * (mu[0] + mut[0] + mu[1] + mut[1]);
    A(0, 0) = mu_lo * (-1. / y[0]) - mu_hi * (1. / (y[1] - y[0]));
    A(0, 1) = mu_hi * (1. / (y[1] - y[0]));
    b(0) = -mu_lo * uw / y[0];
    // freestream BC
    A(n_grid - 1, n_grid - 1) = 1.0;
    b(n_grid - 1) = uwm;
    // interior points
    for (int i = 1; i < n_grid - 1; ++i) {
      mu_lo = 0.5 * (mu[i - 1] + mut[i - 1] + mu[i] + mut[i]);
      mu_hi = 0.5 * (mu[i] + mut[i] + mu[i + 1] + mut[i + 1]);
      A(i, i - 1) = mu_lo / (y[i] - y[i - 1]);
      A(i, i) = -mu_lo / (y[i] - y[i - 1]) - mu_hi / (y[i + 1] - y[i]);
      A(i, i + 1) = mu_hi / (y[i + 1] - y[i]);
      b(i) = 0.0;
    }
    // printMatrix<n_grid>(A, "A");
    // printArray<n_grid>(b, "b");

    // solve
    lusolver.define(A);
    lusolver(tmp, b.arr);

    res_u = 0.0;
    for (int i = 0; i < n_grid; ++i) {
      res_u += std::abs(u[i] - tmp[i]) / Real(n_grid);
      u[i] = tmp[i];
    }
    // printArray(u, n_grid, "u");

    // solve for T
    for (int i = 0; i < n_grid; ++i) {
      for (int j = 0; j < n_grid; ++j) {
        A(i, j) = 0.0;
      }
      b(i) = 0.0;
    }
    // wall BC
    mu_lo = mu[0] + mut[0];
    mu_hi = 0.5 * (mu[0] + mut[0] + mu[1] + mut[1]);
    Real lam_lo = (lam[0] + cp[0] * mut[0] / Prt);
    Real lam_hi =
      0.5 * (lam[0] + cp[0] * mut[0] / Prt + lam[1] + cp[1] * mut[1] / Prt);
    if (T_wall > 0.0) {
      // Isothermal wall
      A(0, 0) = lam_lo * (-1. / y[0]) - lam_hi * (1. / (y[1] - y[0]));
      A(0, 1) = lam_hi * (1. / (y[1] - y[0]));
      b(0) = -lam_lo * T_wall / y[0] + mu_lo * (uw * u[0] / y[0]) -
             mu_hi * ((u[0] + u[1]) / 2 * (u[1] - u[0]) / (y[1] - y[0]));
    } else {
      // adiabatic wall
      A(0, 0) = -lam_hi * (1.0 / (y[1] - y[0]));
      A(0, 1) = lam_hi * (1. / (y[1] - y[0]));
      b(0) = mu_lo * (uw * u[0] / y[0]) -
             mu_hi * ((u[0] + u[1]) / 2 * (u[1] - u[0]) / (y[1] - y[0]));
    }
    // freestream BC
    A(n_grid - 1, n_grid - 1) = 1.0;
    b(n_grid - 1) = Twm;
    // interior points
    for (int i = 1; i < n_grid - 1; ++i) {
      mu_lo = 0.5 * (mu[i - 1] + mut[i - 1] + mu[i] + mut[i]);
      mu_hi = 0.5 * (mu[i] + mut[i] + mu[i + 1] + mut[i + 1]);
      lam_lo = 0.5 * (lam[i - 1] + cp[i - 1] * mut[i - 1] / Prt + lam[i] +
                      cp[i] * mut[i] / Prt);
      lam_hi = 0.5 * (lam[i] + cp[i] * mut[i] / Prt + lam[i + 1] +
                      cp[i + 1] * mut[i + 1] / Prt);
      A(i, i - 1) = lam_lo / (y[i] - y[i - 1]);
      A(i, i) = -lam_lo / (y[i] - y[i - 1]) - lam_hi / (y[i + 1] - y[i]);
      A(i, i + 1) = lam_hi / (y[i + 1] - y[i]);
      b(i) =
        mu_lo * ((u[i - 1] + u[i]) / 2 * (u[i] - u[i - 1]) / (y[i] - y[i - 1])) -
        mu_hi * ((u[i] + u[i + 1]) / 2 * (u[i + 1] - u[i]) / (y[i + 1] - y[i]));
    }
    // printMatrix<n_grid>(A, "C");
    // printArray<n_grid>(b, "d");

    // solve
    lusolver.define(A);
    lusolver(tmp, b.arr);
    res_T = 0.0;
    for (int i = 0; i < n_grid; ++i) {
      res_T += std::abs(T[i] - tmp[i]) / Real(n_grid);
      T[i] = std::max(tmp[i], 90.0);
    }
    // printArray(T, n_grid, "T");

    // std::cout << iter << " res_u = " << res_u << ", res_T = " << res_T << std::endl;
    iter++;
  }

  Real mu0, lam0, dummy;
  if (T_wall > 0.0) {
    trans.transport(false, true, true, false, false, T_wall, rhowm * Twm / T_wall, Y,
                    nullptr, nullptr, mu0, dummy, lam0, ltransparm);
    tau = mu0 * (u[0] - uw) / y[0];
    q = lam0 * (T[0] - T_wall) / y[0];
  } else {
    trans.transport(false, true, false, false, false, T[0], rhowm * Twm / T[0], Y,
                    nullptr, nullptr, mu0, dummy, lam0, ltransparm);
    tau = mu0 * (u[0] - uw) / y[0];
    q = 0.0;
  }

#ifdef WALLMODEL_DEBUG
  trans_parms.deallocate();
#endif
}