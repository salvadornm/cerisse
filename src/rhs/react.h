#ifndef REACT_H
#define REACT_H

#include <PelePhysics.H>
#include <ReactorBase.H>

template <typename cls_t>
class reactor_t {
  std::unique_ptr<pele::physics::reactions::ReactorBase> m_reactor;

 public:
  reactor_t() {
    // These needs amrex::ParmParse, but they are called before main, so before amrex::Initialize
    // std::string reactor_type = "ReactorNull";
    // {
    //   amrex::ParmParse pp("cns");
    //   pp.query("reactor_type", reactor_type);
    // }
    // m_reactor = pele::physics::reactions::ReactorBase::create(reactor_type);
    // m_reactor->init(1, 1);
  }

  ~reactor_t() { 
    // m_reactor->close(); 
  }

  /**
   * @brief Calculate chemical reaction source term, adding to the right-hand-side
   * (rhs) array.
   *
   * @tparam cls_t The problem closure class typename.
   * @param mfi    The MFIter object representing the current grid patch.
   * @param prims  The input primitive variables array.
   * @param rhs    The output array where the updated right-hand-side will be stored.
   * @param cls    The problem closure object (for indicies).
   * @param dt     The time step size.
   */
  void inline src(const amrex::MFIter& mfi,
                  const amrex::Array4<const amrex::Real>& prims,
                  const amrex::Array4<amrex::Real>& rhs, const cls_t& cls,
                  const amrex::Real& dt);
};

// see https://www.codeproject.com/Articles/48575/How-to-Define-a-Template-Class-in-a-h-File-and-Imp
template <typename cls_t>
void inline reactor_t<cls_t>::src(const amrex::MFIter& mfi,
                                  const amrex::Array4<const amrex::Real>& prims,
                                  const amrex::Array4<amrex::Real>& rhs,
                                  const cls_t& cls, const amrex::Real& dt) {
  BL_PROFILE("reactor_t::src()");

  using amrex::Real;
  using amrex::Box;
  using amrex::Array4;
  using amrex::FArrayBox;
  using amrex::IArrayBox;

  const Box bx = mfi.tilebox();

  // TODO: do not work in fine covered box
  // TODO: stochastic fields indexing

  ///////////////////// Prepare for react /////////////////////
  FArrayBox tempf(bx, NUM_SPECIES + 4, The_Async_Arena());
  auto const& rY = tempf.array(0);
  auto const& rEi = tempf.array(NUM_SPECIES);
  auto const& T = tempf.array(NUM_SPECIES + 1);
  auto const& rYsrc = Array4<Real>(rhs, cls.UFS, NUM_SPECIES);
  auto const& rEisrc = tempf.array(NUM_SPECIES + 2);
  auto const& fc =
      tempf.array(NUM_SPECIES + 3);  // number of RHS eval (not used)
  IArrayBox maskf(bx, 1, The_Async_Arena());
  maskf.setVal(1);
  auto const& mask = maskf.array();  // 1: do reaction, -1: skip reaction

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    // [rY, rEi, T, rEisrc]
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
      rY(i, j, k, ns) = prims(i, j, k, cls.QRHO) * prims(i, j, k, cls.QFS + ns);
    }

    rEi(i, j, k) = prims(i, j, k, cls.QRHO) * prims(i, j, k, cls.QEINT);

    T(i, j, k) = prims(i, j, k, cls.QT);

    rEisrc(i, j, k) = rhs(i, j, k, cls.UET);

    Real rho = prims(i, j, k, cls.QRHO);
    Real mx = rho * prims(i, j, k, cls.QU);
    Real my = rho * prims(i, j, k, cls.QV);
    Real mz = rho * prims(i, j, k, cls.QW);
    Real rke = Real(0.5) * (mx * mx + my * my + mz * mz) / rho;
    rho += rhs(i, j, k, cls.URHO) * dt;
    mx += rhs(i, j, k, cls.UMX) * dt;
    my += rhs(i, j, k, cls.UMY) * dt;
    mz += rhs(i, j, k, cls.UMZ) * dt;
    Real rke_new = Real(0.5) * (mx * mx + my * my + mz * mz) / rho;
    rEisrc(i, j, k) = rhs(i, j, k, cls.UET) - (rke_new - rke) / dt;

    // TODO: fill mask
    // mask(i, j, k) = (T(i, j, k) > min_react_temp) ? 1 : -1;
  });

  /////////////////////////// React ///////////////////////////
  Real current_time = 0.0;
//   m_reactor->react(bx, rY, rYsrc, T, rEi, rEisrc, fc, mask, dt, current_time
// #ifdef AMREX_USE_GPU
//                    , amrex::Gpu::gpuStream()
// #endif
//   );
  amrex::Gpu::Device::streamSynchronize();  // TODO: is this necessary?

  //////////////////////// Unpack data ////////////////////////
  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    if (mask(i, j, k) != -1) {
      // Monitor problem cell, do not add reaction source
      bool any_rY_unbounded = false;
      bool temp_below_zero = T(i, j, k) < 0.0;
      for (int ns = 0; ns < NUM_SPECIES; ++ns) {
        any_rY_unbounded |=
            (rY(i, j, k, ns) < -1e-5 || rY(i, j, k, ns) > 1.0 + 1e-5 ||
             std::isnan(rY(i, j, k, ns)));
      }

      if (any_rY_unbounded || temp_below_zero) {
        std::cout << "Post-reaction rY=[ ";
        for (int n = 0; n < NUM_SPECIES; ++n)
          std::cout << rY(i, j, k, n) << " ";
        std::cout << "], T=" << T(i, j, k) << " @ " << i << "," << j << "," << k
                  << '\n';
      } else {
        // Update species source terms
        for (int ns = 0; ns < NUM_SPECIES; ++ns) {
          Real rY_init =
              prims(i, j, k, cls.QRHO) * prims(i, j, k, cls.QFS + ns);
          rhs(i, j, k, cls.UFS + ns) += (rY(i, j, k, ns) - rY_init) / dt;
        }
        // rE is unchanged because it includes chemical energy
      }
    }
  });

  // TODO: Record runtime for load balancing
}

#endif