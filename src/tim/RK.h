#ifndef RK_H_
#define RK_H_

// instead of RK can we generalise to single step explicit?
template <int order, int stages, typename rhs_t>
class RK
{
public:

  void inline static advance(AmrLevel& amrlevel, MultiFab& s1, MultiFab& s2, MultiFab& stemp, FluxRegister* fr_crse, FluxRegister* fr_fine, Real& time, Real& dt, int& ncons, int& nghost)
  {
    if constexpr (order == 0) { // returns rhs
      // FillPatch(amrlevel, stemp, nghost, time, 0, 0, ncons);
      // rhs_t::compute(stemp, dt, fr_crse, fr_fine);
      rhs_t::compute();
      MultiFab::Copy(s2, stemp, 0, 0, ncons, 0);
    };
  //   } else if (order_rk == 1) {
  //     FillPatch(*this, stemp, nghost, time, State_Type, 0,
  //               ncons); // filled at t_n to evalulate f(t_n,y_n).
  //     compute_rhs(stemp, dt, fr_as_crse, fr_as_fine);
  //     MultiFab::LinComb(s2, Real(1.0), s1, 0, dt, stemp, 0, 0, ncons, 0);
  //   } else if (order_rk == 2) {
  //     // Low storage SSPRKm2 with m stages (C = m-1, Ceff=1-1/m). Where C is the
  //     // SSPRK coefficient, it also represents the max CFL over the whole
  //     // integration step (including m stages). From pg 84 Strong Stability
  //     // Preserving Runge–kutta And Multistep Time Discretizations
  //     int m = stages_rk;
  //     // Copy to s2 from s1
  //     MultiFab::Copy(s2, s1, 0, 0, ncons, 0);
  //     state[0].setOldTimeLevel(time);
  //     state[0].setNewTimeLevel(time);
  //     // first to m-1 stages
  //     // Print() << "-------- before RK stages -------" << std::endl;
  //     // Print() << "time = " << time << std::endl;
  //     // Print() << "dt = " << dt << std::endl;
  //     // state[0].printTimeInterval(std::cout);
  //     // Print() << "---------------------------------" << std::endl;
  //     for (int i = 1; i <= m - 1; i++) {
  //       FillPatch(*this, stemp, nghost, time + dt * Real(i - 1) / (m - 1),
  //                 State_Type, 0, ncons);
  //       compute_rhs(stemp, dt / Real(m - 1), fr_as_crse, fr_as_fine);
  //       MultiFab::Saxpy(s2, dt / Real(m - 1), stemp, 0, 0, ncons, 0);
  //       state[State_Type].setNewTimeLevel(
  //         time +
  //         dt * Real(i) / (m - 1)); // important to do this for correct fillpatch
  //                                  // interpolations for the proceeding stages
  //     }
  //     // final stage
  //     FillPatch(*this, stemp, nghost, time + dt, State_Type, 0, ncons);
  //     compute_rhs(stemp, dt / Real(m - 1), fr_as_crse, fr_as_fine);
  //     MultiFab::LinComb(s2, Real(m - 1), s2, 0, dt, stemp, 0, 0, ncons, 0);
  //     MultiFab::LinComb(s2, Real(1.0) / m, s1, 0, Real(1.0) / m, s2, 0, 0, ncons, 0);

  //     state[State_Type].setNewTimeLevel(
  //       time + dt); // important to do this for correct fillpatch
  //                   // interpolations for the proceeding stages

  //     // Print() << "--------- after RK stages --------" << std::endl;
  //     // Print() << "time = " << time << std::endl;
  //     // Print() << "dt = " << dt << std::endl;
  //     // state[0].printTimeInterval(std::cout);
  //     // Print() << "----------------------------------" << std::endl;
  //   } else if (order_rk == 3) {
  //     if (stages_rk == 3) {
  //       state[0].setOldTimeLevel(time);
  //       // http://ketch.github.io/numipedia/methods/SSPRK33.html
  //       // state[0].setOldTimeLevel (time);
  //       FillPatch(*this, stemp, nghost, time, State_Type, 0,
  //                 ncons); // filled at t_n to evalulate f(t_n,y_n).
  //       compute_rhs(stemp, dt, fr_as_crse, fr_as_fine);
  //       MultiFab::LinComb(s2, Real(1.0), s1, 0, dt, stemp, 0, 0, ncons, 0);

  //       state[0].setNewTimeLevel(
  //         time + dt); // same time as upcoming FillPatch ensures we copy s2 to
  //                     // Sborder, without time interpolation
  //       FillPatch(*this, stemp, nghost, time + dt, State_Type, 0, ncons);
  //       compute_rhs(stemp, dt / 4, fr_as_crse, fr_as_fine);
  //       MultiFab::Xpay(stemp, dt, s2, 0, 0, ncons, 0);
  //       MultiFab::LinComb(s2, Real(3.0) / 4, s1, 0, Real(1.0) / 4, stemp, 0, 0,
  //                         ncons, 0);

  //       state[0].setNewTimeLevel(
  //         time + dt / 2); // same time as upcoming FillPatch ensures we copy s2
  //                         // to Sborder, without time interpolation
  //       FillPatch(*this, stemp, nghost, time + dt / 2, State_Type, 0, ncons);
  //       compute_rhs(stemp, dt * Real(2.0) / 3, fr_as_crse, fr_as_fine);
  //       MultiFab::Xpay(stemp, dt, s2, 0, 0, ncons, 0);
  //       MultiFab::LinComb(s2, Real(1.0) / 3, s1, 0, Real(2.0) / 3, stemp, 0, 0,
  //                         ncons, 0);

  //       state[State_Type].setNewTimeLevel(
  //         time + dt); // important to do this for correct fillpatch
  //                     // interpolations for the proceeding stages
  //     } else if (stages_rk == 4) {
  //       // http://ketch.github.io/numipedia/methods/SSPRK43.html and From pg 85
  //       // Strong Stability Preserving Runge–kutta And Multistep Time
  //       // Discretizations

  //       state[0].setOldTimeLevel(time);
  //       FillPatch(*this, stemp, nghost, time, State_Type, 0, ncons);
  //       compute_rhs(stemp, dt / 2, fr_as_crse, fr_as_fine);
  //       MultiFab::LinComb(s2, Real(1.0), s1, 0, dt / 2, stemp, 0, 0, ncons, 0);

  //       state[0].setNewTimeLevel(
  //         time + dt / 2); // same time as upcoming FillPatch ensures we copy s2
  //                         // to Sborder, without time interpolation
  //       FillPatch(*this, stemp, nghost, time + dt / 2, State_Type, 0, ncons);
  //       compute_rhs(stemp, dt / 2, fr_as_crse, fr_as_fine);
  //       MultiFab::Saxpy(s2, dt / 2, stemp, 0, 0, ncons, 0);

  //       state[0].setNewTimeLevel(
  //         time + dt); // same time as upcoming FillPatch ensures we copy s2 to
  //                     // Sborder, without time interpolation
  //       FillPatch(*this, stemp, nghost, time + dt, State_Type, 0, ncons);
  //       compute_rhs(stemp, dt / 6, fr_as_crse, fr_as_fine);
  //       MultiFab::LinComb(s2, Real(2.0) / 3, s1, 0, Real(1.0) / 3, s2, 0, 0, ncons,
  //                         0);
  //       MultiFab::Saxpy(s2, dt / 6, stemp, 0, 0, ncons, 0);

  //       state[0].setNewTimeLevel(
  //         time + dt / 2); // same time as upcoming FillPatch ensures we copy s2
  //                         // to Sborder, without time interpolation
  //       FillPatch(*this, stemp, nghost, time + dt / 2, State_Type, 0, ncons);
  //       compute_rhs(stemp, dt / 2, fr_as_crse, fr_as_fine);
  //       MultiFab::Saxpy(s2, dt / 2, stemp, 0, 0, ncons, 0);

  //       state[State_Type].setNewTimeLevel(
  //         time + dt); // important to do this for correct fillpatch
  //                     // interpolations for the proceeding stages
  //     } else {
  //       // Low storage SSPRKm3 with m=n^2, n>=3 stages (C=2, Ceff=0.5). From pg 85
  //       // Strong Stability Preserving Runge–kutta And Multistep Time
  //       // Discretizations
  //       // TODO Generally SSPRK(n^2,3) where n>2 - Ceff=1-1/n
  //       Print() << "SSPRK(m^2)3 not implemented yet" << std::endl;
  //       exit(0);
  //     }
  //   }
  };
};

// template
// class implicit{
// }

#endif