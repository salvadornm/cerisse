#ifndef CNS_LES_H_
#define CNS_LES_H_

#include <Factory.H>

#include <AMReX_REAL.H>
#include <AMReX_Array4.H>

#include "index_macros.H"

using namespace amrex;

/**
 * @brief Base class to register eddy viscosity type LES models.
 */
class LESModel : public pele::physics::Factory<LESModel>
{
public:
  static const std::string base_identifier() { return "LESModel"; }

  ~LESModel() override = default;

  virtual void mu_T_cc (const int i, const int j, const int k,
                        const Array4<const Real>& q,
                        const GpuArray<Real, AMREX_SPACEDIM>& dxinv,
                        const Real deltabar,
                        const Real Cs,
                        Real& mu_T) = 0;
  
  virtual void mu_T_fc (const Real dudx, const Real dudy, const Real dudz, 
                        const Real dvdx, const Real dvdy, const Real dvdz, 
                        const Real dwdx, const Real dwdy, const Real dwdz, 
                        const Real rho,
                        const Real deltabar,
                        const Real Cs,
                        Real& mu_T) = 0;  
};

//
// Declare LES models

class Smagorinsky : public LESModel::Register<Smagorinsky>
{
public:
  static const std::string identifier() { return "Smagorinsky"; }

  void mu_T_cc (const int i, const int j, const int k,
                const Array4<const Real>& q,
                const GpuArray<Real, AMREX_SPACEDIM>& dxinv,
                const Real deltabar,
                const Real Cs,
                Real& mu_T) override;
  
  void mu_T_fc (const Real dudx, const Real dudy, const Real dudz, 
                const Real dvdx, const Real dvdy, const Real dvdz, 
                const Real dwdx, const Real dwdy, const Real dwdz, 
                const Real rho,
                const Real deltabar,
                const Real Cs,
                Real& mu_T) override;
};

// Nicoud & Ducros (1999) https://doi.org/10.1023/A:1009995426001
class WALE : public LESModel::Register<WALE>
{
public:
  static const std::string identifier() { return "WALE"; }

  void mu_T_cc (const int i, const int j, const int k,
                const Array4<const Real>& q,
                const GpuArray<Real, AMREX_SPACEDIM>& dxinv,
                const Real deltabar,
                const Real Cs,
                Real& mu_T) override;
  
  void mu_T_fc (const Real dudx, const Real dudy, const Real dudz, 
                const Real dvdx, const Real dvdy, const Real dvdz, 
                const Real dwdx, const Real dwdy, const Real dwdz, 
                const Real rho,
                const Real deltabar,
                const Real Cs,
                Real& mu_T) override;
};

#endif