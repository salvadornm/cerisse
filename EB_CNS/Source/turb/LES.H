#ifndef CNS_LES_H_
#define CNS_LES_H_

#include <AMReX_Array4.H>
#include <AMReX_REAL.H>

#include "index_macros.H"

/**
 * @brief Base class to register eddy viscosity type LES models.
 */
class LESModel
{
public:
  static const std::string base_identifier() { return "LESModel"; }

  static void create(std::string selector, LESModel*& model);
  
  static void destroy(LESModel*& model);

  virtual AMREX_GPU_DEVICE void mu_T_cc(
    const int i, const int j, const int k, const amrex::Array4<const amrex::Real>& q,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxinv,
    const amrex::Real deltabar, const amrex::Real Cs, amrex::Real& mu_T) = 0;

  virtual AMREX_GPU_DEVICE void
  mu_T_fc(const amrex::Real dudx, const amrex::Real dudy, const amrex::Real dudz,
          const amrex::Real dvdx, const amrex::Real dvdy, const amrex::Real dvdz,
          const amrex::Real dwdx, const amrex::Real dwdy, const amrex::Real dwdz,
          const amrex::Real rho, const amrex::Real deltabar, const amrex::Real Cs,
          amrex::Real& mu_T) = 0;

  virtual AMREX_GPU_DEVICE void xi_T_cc(
    const int i, const int j, const int k, const amrex::Array4<const amrex::Real>& q,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxinv,
    const amrex::Real deltabar, const amrex::Real CI, amrex::Real& xi_T);
};

//
// Declare LES models

class Smagorinsky : public LESModel
{
public:
  static const std::string identifier() { return "Smagorinsky"; }

  AMREX_GPU_DEVICE void mu_T_cc(
    const int i, const int j, const int k, const amrex::Array4<const amrex::Real>& q,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxinv,
    const amrex::Real deltabar, const amrex::Real Cs, amrex::Real& mu_T) override;

  AMREX_GPU_DEVICE void
  mu_T_fc(const amrex::Real dudx, const amrex::Real dudy, const amrex::Real dudz,
          const amrex::Real dvdx, const amrex::Real dvdy, const amrex::Real dvdz,
          const amrex::Real dwdx, const amrex::Real dwdy, const amrex::Real dwdz,
          const amrex::Real rho, const amrex::Real deltabar, const amrex::Real Cs,
          amrex::Real& mu_T) override;
};

// Nicoud & Ducros (1999) https://doi.org/10.1023/A:1009995426001
class WALE : public LESModel
{
public:
  static const std::string identifier() { return "WALE"; }

  AMREX_GPU_DEVICE void mu_T_cc(
    const int i, const int j, const int k, const amrex::Array4<const amrex::Real>& q,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxinv,
    const amrex::Real deltabar, const amrex::Real Cs, amrex::Real& mu_T) override;

  AMREX_GPU_DEVICE void
  mu_T_fc(const amrex::Real dudx, const amrex::Real dudy, const amrex::Real dudz,
          const amrex::Real dvdx, const amrex::Real dvdy, const amrex::Real dvdz,
          const amrex::Real dwdx, const amrex::Real dwdy, const amrex::Real dwdz,
          const amrex::Real rho, const amrex::Real deltabar, const amrex::Real Cs,
          amrex::Real& mu_T) override;
};

#endif