#ifndef CNS_WALL_MODEL_H_
#define CNS_WALL_MODEL_H_

#include <AMReX_REAL.H>
#include <Factory.H>
#include <PelePhysics.H>

#include "index_macros.H"

class LESWallModel : public pele::physics::Factory<LESWallModel>
{
public:
  static const std::string base_identifier() { return "LESWallModel"; }

  ~LESWallModel() override = default;

  /**
   * @brief Compute the parallel wall stress and heat flux
   *
   * @param u Parallel velocity at coupling point (y = h)
   * @param T Temperature at coupling point (y = h)
   * @param p Pressure (assumed to be constant along wall normal direction)
   * @param Y Species mass fractions (assmed to be constant as well)
   * @param h Coupling point height
   * @param T_wall Temperature at the wall (set to <0 if adiaiabatic)
   * @param tau[out] Parallel wall stress
   * @param q[out]   Wall heat flux
   */
  virtual void parallel_wall_stress(amrex::Real u, amrex::Real T, amrex::Real rho,
                                    amrex::Real Y[NUM_SPECIES], amrex::Real h,
                                    amrex::Real mu, amrex::Real lam,
                                    amrex::Real T_wall, amrex::Real& tau,
                                    amrex::Real& q) = 0;
};

// Algebric wall model

class LawOfTheWall : public LESWallModel::Register<LawOfTheWall>
{
public:
  static const std::string identifier() { return "LawOfTheWall"; }

  void parallel_wall_stress(amrex::Real u, amrex::Real T, amrex::Real rho,
                            amrex::Real Y[NUM_SPECIES], amrex::Real h,
                            amrex::Real mu, amrex::Real lam, amrex::Real T_wall,
                            amrex::Real& tau, amrex::Real& q) override;
};

// ODE wall model

class EquilibriumODE : public LESWallModel::Register<EquilibriumODE>
{
public:
  static const std::string identifier() { return "EquilibriumODE"; }

  void parallel_wall_stress(amrex::Real u, amrex::Real T, amrex::Real rho,
                            amrex::Real Y[NUM_SPECIES], amrex::Real h,
                            amrex::Real mu, amrex::Real lam, amrex::Real T_wall,
                            amrex::Real& tau, amrex::Real& q) override;
};

// Apply wall model on regular walls

template <class WallModel>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void apply_regular_wall_model(
  const amrex::IntVect& iv, const int dir, const bool lo_is_wall,
  const bool hi_is_wall, const int domlo, const int domhi,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx,
  const amrex::Array4<const amrex::Real>& q,
  const amrex::Array4<const amrex::Real>& coefs, const bool is_isothermal,
  const amrex::Real wall_temp, amrex::Real flx[UFS])
{
  if ((iv[dir] == domhi + 1 && hi_is_wall) || (iv[dir] == domlo && lo_is_wall)) {
    // first fluid cell (1)
    const auto ivd = amrex::IntVect::TheDimensionVector(dir);
    const auto iv1 = (iv[dir] == domlo) ? iv : iv - ivd;

    // get density and transport coefficients
    // assume adiabatic wall (this is just to estimate y+)
    amrex::Real Yw[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; ++n) { Yw[n] = q(iv1, QFS + n); }
    const amrex::Real rhow = q(iv1, QRHO);
    const amrex::Real muw = coefs(iv1, CMU);
    const amrex::Real lamw = coefs(iv1, CLAM);

    // check y+
    const amrex::Real local_dx = dx[dir];
    const amrex::Real tau =
      std::sqrt(flx[UMX] * flx[UMX] + flx[UMY] * flx[UMY] + flx[UMZ] * flx[UMZ]);
    const amrex::Real utau = std::sqrt(tau / rhow);
    const amrex::Real yplus = rhow * utau / muw * local_dx;

    // setup image point (2) values
    // image point (second fluid point)
    const auto iv2 = (iv[dir] == domlo) ? iv1 + ivd : iv1 - ivd;
    const amrex::Real T2 = q(iv2, QTEMP);
    const amrex::Real rho2 = q(iv2, QRHO);
    // tangent vectors (t1.n = 0, t2.t1 = 0, u.t2 = 0)
    amrex::Real u2[3] = {q(iv2, QU), q(iv2, QV), q(iv2, QW)};
    u2[dir] = 0.0; // wall-normal component = 0
    const amrex::Real u2_parallel =
      std::sqrt(u2[0] * u2[0] + u2[1] * u2[1] + u2[2] * u2[2]);
    const amrex::Real t1[3] = {u2[0] / u2_parallel, u2[1] / u2_parallel,
                               u2[2] / u2_parallel};

    // call wm.parallel_wall_stress() and overwrite wall fluxes
    // prevent numerical issues at u_parallel ~ 0 and extreme temperatures
    if (yplus > 11.0 && u2_parallel > 100.0 && !(T2 < 90.0) && !(T2 > 4000.0) &&
        !(rho2 < 0.0)) {
      const amrex::Real T_wall = is_isothermal ? wall_temp : -1.0;
      const amrex::Real h = 1.5 * local_dx;
      amrex::Real tauw, qw;
      WallModel wm;
      wm.parallel_wall_stress(u2_parallel, T2, rho2, Yw, h, muw, lamw, T_wall, tauw,
                              qw);
      if (!std::isnan(tauw) && !std::isnan(qw)) {
        // this accouts for the opposite normal vector directions in lo and hi walls
        const amrex::Real lohi = (iv[dir] == domlo) ? -1.0 : 1.0;
        flx[UMX] = t1[0] * tauw * lohi;
        flx[UMY] = t1[1] * tauw * lohi;
        flx[UMZ] = t1[2] * tauw * lohi;
        flx[UEDEN] = qw * lohi;
      }
    }
  }
}

template <class WallModel>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void apply_regular_wall_model_sf(
  const amrex::IntVect& iv, const int dir, const bool lo_is_wall,
  const bool hi_is_wall, const int domlo, const int domhi,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx,
  const amrex::Array4<const amrex::Real>& q0, const bool is_isothermal,
  const amrex::Real wall_temp,
  const pele::physics::transport::TransParm<
    pele::physics::EosType, pele::physics::PhysicsType::transport_type>* ltransparm,
  const amrex::Array4<amrex::Real>& flx)
{
  if ((iv[dir] == domhi + 1 && hi_is_wall) || (iv[dir] == domlo && lo_is_wall)) {
    // first fluid cell (1)
    const auto ivd = amrex::IntVect::TheDimensionVector(dir);
    const auto iv1 = (iv[dir] == domlo) ? iv : iv - ivd;

    // get density and transport coefficients
    amrex::Real Yw[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; ++n) { Yw[n] = q0(iv1, 5 + n); }
    const amrex::Real rhow = q0(iv1, 0);
    auto trans = pele::physics::PhysicsType::transport();
    // auto const* ltransparm = CNS::trans_parms.device_trans_parm();
    amrex::Real muw, lamw, xi_unused;
    trans.transport(true, true, true, false, false, wall_temp, rhow, Yw, nullptr,
                    nullptr, muw, xi_unused, lamw, ltransparm);

    // check y+
    const amrex::Real local_dx = dx[dir];
    const amrex::Real tau =
      std::sqrt(flx(iv, UMX) * flx(iv, UMX) + flx(iv, UMY) * flx(iv, UMY) +
                flx(iv, UMZ) * flx(iv, UMZ)) /
      amrex::Real(NUM_FIELD);
    const amrex::Real yplus = std::sqrt(tau * rhow) / muw * local_dx;

    // setup image point (2) values
    // image point (second fluid point)
    const auto iv2 = iv[dir] == domlo ? iv1 + ivd : iv1 - ivd;
    const amrex::Real T2 = q0(iv2, 4);
    const amrex::Real rho2 = q0(iv2, 0);
    // tangent vectors (t1.n = 0, t2.t1 = 0, u.t2 = 0)
    amrex::Real u2[3] = {q0(iv2, 1), q0(iv2, 2), q0(iv2, 3)};
    u2[dir] = 0.0;
    const amrex::Real u_parallel =
      std::sqrt(u2[0] * u2[0] + u2[1] * u2[1] + u2[2] * u2[2]);
    const amrex::Real t1[3] = {u2[0] / u_parallel, u2[1] / u_parallel,
                               u2[2] / u_parallel};

    // call wm.parallel_wall_stress() and overwrite wall fluxes
    // prevent numerical issues at u_parallel ~ 0 and extreme temperatures
    if (yplus > 11.0 && u_parallel > 100.0 && !(T2 < 90.0) && !(T2 > 4000.0) &&
        !(rho2 < 0.0)) {
      WallModel wm;
      amrex::Real T_wall = is_isothermal ? wall_temp : -1.0;
      amrex::Real h = 1.5 * local_dx;
      amrex::Real tauw, qw;
      wm.parallel_wall_stress(u_parallel, T2, rho2, Yw, h, muw, lamw, T_wall, tauw,
                              qw);
      if (!std::isnan(tauw) && !std::isnan(qw)) {
        const amrex::Real lohi = iv[dir] == domlo ? -1.0 : 1.0;
        for (int nf = 1; nf <= NUM_FIELD; ++nf) {
          flx(iv, UMX) = t1[0] * tauw * lohi * amrex::Real(NUM_FIELD);
          flx(iv, UMY) = t1[1] * tauw * lohi * amrex::Real(NUM_FIELD);
          flx(iv, UMZ) = t1[2] * tauw * lohi * amrex::Real(NUM_FIELD);
          flx(iv, UEDEN) = qw * lohi * amrex::Real(NUM_FIELD);
        }
      }
    }
  }
}

#endif