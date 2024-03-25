#ifndef REACT_H
#define REACT_H

#include <PelePhysics.H>
#include <ReactorBase.H>

// template <int reactor_type, typename cls_t>
template <typename cls_t>
class reactor_t {
 public:
 public:
  bool m_initialized = false;
  std::unique_ptr<pele::physics::reactions::ReactorBase> m_reactor;

  // reactor_t() {
  //   if (reactor_type == 1) {
  //     m_reactor = pele::physics::reactions::ReactorBase::create("ReactorCvode");
  //     m_initialized = true;
  //   } else {
  //     amrex::Abort("reactor_t not initialised");
  //   }
  //   m_reactor->init(1, 1);
  // };

  reactor_t() {
    std::string reactor_type;
    {
      amrex::ParmParse pp("cns");
      pp.get("reactor_type", reactor_type);
    }
    m_reactor = pele::physics::reactions::ReactorBase::create(reactor_type);
    m_reactor->init(1, 1);
  };

  ~reactor_t() {
    if (m_initialized) m_reactor->close();
  }

  /**
   * @brief Calculate chemical reaction source term, adding to the
   * right-hand-side (rhs) array.
   *
   * @tparam cls_t The problem closure class typename.
   * @param mfi    The MFIter object representing the current grid patch.
   * @param prims  The input primitive variables array.
   * @param rhs    The output array where the updated right-hand-side will be
   * stored.
   * @param cls    The problem closure object (for indicies).
   * @param dt     The time step size. (react() requires it to be non-const)
   */
  // https://www.codeproject.com/Articles/48575/How-to-Define-a-Template-Class-in-a-h-File-and-Imp
  void inline src(const amrex::MFIter& mfi,
                  const amrex::Array4<const amrex::Real>& prims,
                  const amrex::Array4<amrex::Real>& rhs, const cls_t* cls_d,
                  amrex::Real dt) {
    if (!m_initialized) amrex::Abort("reactor_t not initialised");

    // amrex::Print() << "reactor_t::src()" << std::endl;

    BL_PROFILE("reactor_t::src()");

    // put here because this is a .h file
    using amrex::Array4;
    using amrex::Box;
    using amrex::FArrayBox;
    using amrex::IArrayBox;
    using amrex::Real;

    const Box bx = mfi.tilebox();

      // TODO: do not work in fine covered box
      // TODO: stochastic fields indexing

      ///////////////////// Prepare for react /////////////////////
      FArrayBox tempf(bx, 2 * NUM_SPECIES + 4, The_Async_Arena());
      auto const& rY = tempf.array(0);
      auto const& rEi = tempf.array(NUM_SPECIES);
      auto const& T = tempf.array(NUM_SPECIES + 1);
      auto const& rYsrc = tempf.array(NUM_SPECIES + 2);
      auto const& rEisrc = tempf.array(2 * NUM_SPECIES + 2);
      auto const& fc =
       
        tempf.array(2 * NUM_SPECIES + 3);  // number of RHS eval (not used)
      IArrayBox maskf(bx, 1, The_Async_Arena());
      maskf.setVal<RunOn::Gpu>(1);
      auto const& mask = maskf.array();  // 1: do reaction, -1: skip reaction

      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const auto& cls = *cls_d;

        // [rY, rEi, T, rYsrc, rEisrc] Remember to convert to CGS!!
        for (int ns = 0; ns < NUM_SPECIES; ++ns) {
          rY(i, j, k, ns) =
            prims(i, j, k, cls.QRHO) * prims(i, j, k, cls.QFS + ns) * 1.0e-3;
          rYsrc(i, j, k, ns) = rhs(i, j, k, cls.UFS + ns) * 1.0e-3;
        }

      rEi(i, j, k) =
          prims(i, j, k, cls.QRHO) * prims(i, j, k, cls.QEINT) * 10.0;
      rEi(i, j, k) =
          prims(i, j, k, cls.QRHO) * prims(i, j, k, cls.QEINT) * 10.0;

      T(i, j, k) = prims(i, j, k, cls.QT);
      AMREX_ALWAYS_ASSERT(T(i, j, k) > 0.0);
      T(i, j, k) = prims(i, j, k, cls.QT);
      AMREX_ALWAYS_ASSERT(T(i, j, k) > 0.0);

      Real rho = prims(i, j, k, cls.QRHO);
      Real mx = rho * prims(i, j, k, cls.QU);
      Real my = rho * prims(i, j, k, cls.QV);
      Real mz = rho * prims(i, j, k, cls.QW);
      Real rke = Real(0.5) * (mx * mx + my * my + mz * mz) / rho;
      for (int ns = 0; ns < NUM_SPECIES; ++ns) {
        rho += rhs(i, j, k, cls.UFS + ns) * dt;
      }
      mx += rhs(i, j, k, cls.UMX) * dt;
      my += rhs(i, j, k, cls.UMY) * dt;
      mz += rhs(i, j, k, cls.UMZ) * dt;
      Real rke_new = Real(0.5) * (mx * mx + my * my + mz * mz) / rho;
      rEisrc(i, j, k) = (rhs(i, j, k, cls.UET) - (rke_new - rke) / dt) * 10.0;
      Real rho = prims(i, j, k, cls.QRHO);
      Real mx = rho * prims(i, j, k, cls.QU);
      Real my = rho * prims(i, j, k, cls.QV);
      Real mz = rho * prims(i, j, k, cls.QW);
      Real rke = Real(0.5) * (mx * mx + my * my + mz * mz) / rho;
      for (int ns = 0; ns < NUM_SPECIES; ++ns) {
        rho += rhs(i, j, k, cls.UFS + ns) * dt;
      }
      mx += rhs(i, j, k, cls.UMX) * dt;
      my += rhs(i, j, k, cls.UMY) * dt;
      mz += rhs(i, j, k, cls.UMZ) * dt;
      Real rke_new = Real(0.5) * (mx * mx + my * my + mz * mz) / rho;
      rEisrc(i, j, k) = (rhs(i, j, k, cls.UET) - (rke_new - rke) / dt) * 10.0;

      // TODO: fill mask
      // mask(i, j, k) = (T(i, j, k) > min_react_temp) ? 1 : -1;
    });
      // TODO: fill mask
      // mask(i, j, k) = (T(i, j, k) > min_react_temp) ? 1 : -1;
    });

    /////////////////////////// React ///////////////////////////
    Real current_time = 0.0;

    // Not necessary to start a stream here, however pelePhysics function only takes a stream. Practically, launch and execution overhead determines  efficiency effect -- https://stackoverflow.com/questions/27038162/how-bad-is-it-to-launch-many-small-kernels-in-cuda#:~:text=Launch%20overhead%3A%20The%20overhead%20of,as%20the%20kernel%20in%20question. Seems unlikely this kernel launch cost will outweigh execution costs.
    /////////////////////////// React ///////////////////////////
    Real current_time = 0.0;

    // Not necessary to start a stream here, however pelePhysics function only takes a stream. Practically, launch and execution overhead determines  efficiency effect -- https://stackoverflow.com/questions/27038162/how-bad-is-it-to-launch-many-small-kernels-in-cuda#:~:text=Launch%20overhead%3A%20The%20overhead%20of,as%20the%20kernel%20in%20question. Seems unlikely this kernel launch cost will outweigh execution costs.
#ifdef AMREX_USE_GPU
    m_reactor->react(bx, rY, rYsrc, T, rEi, rEisrc, fc, mask, dt, current_time,
                     amrex::Gpu::gpuStream());
#else
    m_reactor->react(bx, rY, rYsrc, T, rEi, rEisrc, fc, mask, dt, current_time);
    m_reactor->react(bx, rY, rYsrc, T, rEi, rEisrc, fc, mask, dt, current_time,
                     amrex::Gpu::gpuStream());
#else
    m_reactor->react(bx, rY, rYsrc, T, rEi, rEisrc, fc, mask, dt, current_time);
#endif
    amrex::Gpu::Device::streamSynchronize();  // Important
    amrex::Gpu::Device::streamSynchronize();  // Important

    //////////////////////// Unpack data ////////////////////////
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const auto& cls = *cls_d;
      if (mask(i, j, k) != -1) {
        // Monitor problem cell, do not add reaction source
        bool any_rY_unbounded = false;
        bool temp_leq_zero = (T(i, j, k) <= 0.0);
        for (int ns = 0; ns < NUM_SPECIES; ++ns) {
          any_rY_unbounded |=
              (rY(i, j, k, ns) < -1e-5 || rY(i, j, k, ns) > 1.0 + 1e-5 ||
               std::isnan(rY(i, j, k, ns)));
        }

        if (any_rY_unbounded || temp_leq_zero) {
          printf("Post-reaction rY=[ ");
          for (int ns = 0; ns < NUM_SPECIES; ++ns)
            printf(" %f", rY(i, j, k, ns));
          printf("], T=%f  @ (i,j,k) = (%f, %f, %f)", T(i, j, k), i, j, k);
        } else {
          // Update species source terms
          for (int ns = 0; ns < NUM_SPECIES; ++ns) {
            // rY is overwritten by rY + rYsrc * dt + chem_src * dt = rY +
            // new_rhs * dt Remember to convert from CGS back to SI!!
            Real rY_init =
                prims(i, j, k, cls.QRHO) * prims(i, j, k, cls.QFS + ns);
            rhs(i, j, k, cls.UFS + ns) =
                (rY(i, j, k, ns) * 1.0e3 - rY_init) / dt;
          }
          // rE is unchanged because it includes chemical energy
        }
      }
    });
        if (any_rY_unbounded || temp_leq_zero) {
          printf("Post-reaction rY=[ ");
          for (int ns = 0; ns < NUM_SPECIES; ++ns)
            printf(" %f", rY(i, j, k, ns));
          printf("], T=%f  @ (i,j,k) = (%f, %f, %f)", T(i, j, k), i, j, k);
        } else {
          // Update species source terms
          for (int ns = 0; ns < NUM_SPECIES; ++ns) {
            // rY is overwritten by rY + rYsrc * dt + chem_src * dt = rY +
            // new_rhs * dt Remember to convert from CGS back to SI!!
            Real rY_init =
                prims(i, j, k, cls.QRHO) * prims(i, j, k, cls.QFS + ns);
            rhs(i, j, k, cls.UFS + ns) =
                (rY(i, j, k, ns) * 1.0e3 - rY_init) / dt;
          }
          // rE is unchanged because it includes chemical energy
        }
      }
    });

    // TODO: Record runtime for load balancing
  }
};
    // TODO: Record runtime for load balancing
  }
};

#endif