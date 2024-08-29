#ifndef CNS_WALL_MODEL_H_
#define CNS_WALL_MODEL_H_

#include <AMReX_REAL.H>
#include <Factory.H>

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

#endif