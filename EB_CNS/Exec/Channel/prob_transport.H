#ifdef USE_PROB_TRANSPORT
#ifndef PROBTRANSPORT_H
#define PROBTRANSPORT_H

#include "TransportParams.H"
#include "EOS.H"

#define TRANSPORTTYPES_H

namespace pele {
namespace physics {
namespace transport {
struct ProbTransport;
}
using TransportType = transport::ProbTransport;
} // namespace physics
} // namespace pele


namespace pele {
namespace physics {
namespace transport {

template <typename EOSType>
struct TransParm<EOSType, ProbTransport>
{
  amrex::Real viscosity_exponent{0.7};
  amrex::Real reynolds{0.0};
  // amrex::Real prandtl{0.0};
  amrex::Real mach{0.0};
};

template <typename EOSType>
struct InitTransParm<EOSType, ProbTransport>
{
  void operator()(TransParm<EOSType, ProbTransport>* tparm)
  {
    amrex::ParmParse pp("transport");
    pp.query("viscosity_exponent", tparm->viscosity_exponent);
    pp.query("reynolds", tparm->reynolds);
    pp.query("mach", tparm->mach);
  }
};

struct ProbTransport
{
  using transport_type = ProbTransport;

  static std::string identifier() { return "ProbTransport"; }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE static void transport(
    const bool wtr_get_xi,
    const bool wtr_get_mu,
    const bool wtr_get_lam,
    const bool wtr_get_Ddiag,
    const amrex::Real Tloc,
    const amrex::Real /*rholoc*/,
    amrex::Real* /*Yloc*/,
    amrex::Real* Ddiag,
    amrex::Real& mu,
    amrex::Real& xi,
    amrex::Real& lam,
    TransParm<EosType, transport_type> const* tparm)
  { 
    amrex::Real viscosity = pow(Tloc, tparm->viscosity_exponent);

    if (wtr_get_mu) {
      mu = viscosity / tparm->reynolds;
    }

    if (wtr_get_lam) {
      auto eos = pele::physics::PhysicsType::eos();
      lam = viscosity / tparm->reynolds / tparm->mach / tparm->mach / (eos.gamma - 1);
    }
  
    if (wtr_get_xi) {
      xi = 0.0;
    }

    if (wtr_get_Ddiag) {
      for (int i = 0; i < NUM_SPECIES; ++i) {
        Ddiag[i] = 0.0;
      }
    }
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE static void get_transport_coeffs(
    amrex::Box const& bx,
    amrex::Array4<const amrex::Real> const& /*Y_in*/,
    amrex::Array4<const amrex::Real> const& /*T_in*/,
    amrex::Array4<const amrex::Real> const& /*Rho_in*/,
    amrex::Array4<amrex::Real> const& D_out,
    amrex::Array4<amrex::Real> const& mu_out,
    amrex::Array4<amrex::Real> const& xi_out,
    amrex::Array4<amrex::Real> const& lam_out,
    TransParm<EosType, transport_type> const* tparm)
  {
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    const bool wtr_get_xi = true;
    const bool wtr_get_mu = true;
    const bool wtr_get_lam = true;
    const bool wtr_get_Ddiag = true;

    for (int k = lo.z; k <= hi.z; ++k) {
      for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {

          amrex::Real muloc, xiloc, lamloc;
          amrex::Real Ddiag[NUM_SPECIES] = {0.0};
          transport(
            wtr_get_xi, wtr_get_mu, wtr_get_lam, wtr_get_Ddiag, 0.0, 0.0,
            nullptr, Ddiag, muloc, xiloc, lamloc, tparm);

          // mu, xi and lambda are stored after D in the diffusion multifab
          for (int n = 0; n < NUM_SPECIES; ++n) {
            D_out(i, j, k, n) = Ddiag[n];
          }

          mu_out(i, j, k) = muloc;
          xi_out(i, j, k) = xiloc;
          lam_out(i, j, k) = lamloc;
        }
      }
    }
  }

  template <class... Args>
  AMREX_GPU_HOST_DEVICE explicit ConstTransport(Args... /*unused*/)
  {
  }
};
} // namespace transport
} // namespace physics
} // namespace pele
#endif
#endif