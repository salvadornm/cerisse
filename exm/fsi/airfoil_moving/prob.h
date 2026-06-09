#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>

#include <Closures.h>
#include <RHS.h>

#include <ibm_solver.h>
#include <Constants.h>
#include <ibm_walltypes.h>

using namespace amrex;
using namespace universal_constants;

namespace PROB {

// ============================================================================
// Rigid-body motion: 2D → 3 DOF (tx, ty, θ), 3D → 6 DOF (tx,ty,tz, θx,θy,θz)
// ============================================================================

namespace Motion {
  // Runtime state — updated every timestep
  AMREX_GPU_MANAGED inline Real sim_time = 0.0_rt;

  // Translation oscillation: x(t) = Ax*sin(2π*fx*t), y(t) = Ay*sin(2π*fy*t + phase_y)
  AMREX_GPU_MANAGED inline Real amp_x  = 0.0_rt;   // [m] x-amplitude
  AMREX_GPU_MANAGED inline Real amp_y  = 0.0_rt;   // [m] y-amplitude
  AMREX_GPU_MANAGED inline Real freq_x = 0.0_rt;   // [Hz] x-frequency
  AMREX_GPU_MANAGED inline Real freq_y = 0.0_rt;   // [Hz] y-frequency
  AMREX_GPU_MANAGED inline Real phase_y = 0.0_rt;  // [rad] y-phase offset (π/2 → Lissajous)

  // Rotation oscillation: θ(t) = θ0 + Aθ*sin(2π*fθ*t)
  AMREX_GPU_MANAGED inline Real amp_theta  = 0.0_rt; // [deg] pitch amplitude
  AMREX_GPU_MANAGED inline Real freq_theta = 0.0_rt; // [Hz] pitch frequency

  // Rotation center
  AMREX_GPU_MANAGED inline Real pitch_center_x = 0.0_rt;
  AMREX_GPU_MANAGED inline Real pitch_center_y = 0.0_rt;

  // True motion flag
  AMREX_GPU_MANAGED inline bool true_motion = false;

  static constexpr Real PI = 3.14159265358979323846_rt;

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  Real deg2rad(Real x) { return x * (PI / 180.0_rt); }
}

// ============================================================================
// rigid_trajectory — prescribed 2D rigid-body trajectory (3 DOF)
//
// Returns the displacement (dx, dy) and rotation angle (theta) at a given time.
// All quantities are ABSOLUTE (relative to t=0 reference configuration).
//
// Users can modify this function to implement arbitrary prescribed motion.
// ============================================================================
struct RigidState {
  Real dx;      // x-displacement from reference
  Real dy;      // y-displacement from reference
  Real theta;   // rotation angle [rad] from reference (positive = CCW)
  Real vx;      // x-velocity (for wall BC)
  Real vy;      // y-velocity (for wall BC)
  Real omega;   // angular velocity [rad/s] (for wall BC)
};

AMREX_FORCE_INLINE
RigidState rigid_trajectory(Real time)
{
  using namespace Motion;
  constexpr Real twopi = 2.0_rt * PI;
  RigidState s{};

  // --- Translation: sinusoidal oscillation ---
  //   x(t) = Ax * sin(2π fx t)
  //   y(t) = Ay * sin(2π fy t + phase_y)
  // When phase_y = π/2 and fx = fy → circular orbit (Lissajous)
  s.dx = amp_x * std::sin(twopi * freq_x * time);
  s.dy = amp_y * std::sin(twopi * freq_y * time + phase_y);
  s.vx = amp_x * twopi * freq_x * std::cos(twopi * freq_x * time);
  s.vy = amp_y * twopi * freq_y * std::cos(twopi * freq_y * time + phase_y);

  // --- Rotation: sinusoidal pitch ---
  //   θ(t) = Aθ * sin(2π fθ t)   [deg → rad]
  Real theta_deg = amp_theta * std::sin(twopi * freq_theta * time);
  Real omega_deg = amp_theta * twopi * freq_theta * std::cos(twopi * freq_theta * time);
  s.theta = deg2rad(theta_deg);
  s.omega = deg2rad(omega_deg);

  return s;
}

// constants
static constexpr Real Mw       = 28.96e-3;
static constexpr Real gam      = 1.4;
static constexpr Real Rgas     = gas_constant/Mw;
static constexpr Real Cv       = Rgas/(gam - 1.0);
static constexpr Real Cp       = gam*Cv;

static constexpr bool srp_on   = false;
static constexpr int ibm_eorder= 1;

//////////////////////////// Physical modelling ////////////////////////////////
struct ProbParm
{
  // Quiescent atmosphere — FSI motion test (no freestream)
  static constexpr Real Ma_oo   = 0.0;
  static constexpr Real p_oo    = 101325.0;       // [Pa] sea level
  static constexpr Real T_oo    = 300.0;           // [K]

  // Isothermal wall
  static constexpr Real Twall   = 300.0;

  // Derived
  static constexpr Real rho_oo  = p_oo/(Rgas*T_oo);
  static constexpr Real c_oo    = 347.2;           // sqrt(gamma*Rgas*T_oo)
  static constexpr Real u_oo    = c_oo*Ma_oo;      // = 0
  static constexpr Real eint_oo = Cv*T_oo;
  static constexpr Real kin_oo  = 0.5*rho_oo*u_oo*u_oo;  // = 0

  // Airfoil geometry
  static constexpr Real chord   = 0.1;
  static constexpr Real x_le    = -0.05;
  static constexpr Real x_te    = 0.05;
};

struct methodparm_t {
  public:
  static constexpr int  order = 2;
  static constexpr bool use_LES = false;
};

struct ibmparm_t {
  public:
  static constexpr int  interp_order = 1;
  static constexpr int  extrap_order = ibm_eorder;
  static constexpr int  interp_order_surf = 1;
  static constexpr int  extrap_order_surf = 1;
  static constexpr Real alpha        = 0.4;
  static constexpr Real alpha_surf   = 0.4;
  static constexpr int  ghost_layers = 1;
  static constexpr bool interior_is_solid = true;
};

// CLOSURES
using ProbClosures = closures_dt< indicies_t, transport_suth_t, calorifically_perfect_gas_t<indicies_t> >;

// RHS: WENO-Z5 + viscous
typedef rhs_dt<weno_t<ReconScheme::WenoZ5, ProbClosures>, viscous_t<methodparm_t,ProbClosures>, no_source_t > ProbRHS;

template < typename param, typename cls_t > class ibm_user_t;

typedef ibm_user_t<ProbParm,ProbClosures> TypeWall;
typedef ibm_solver_t<TypeWall,ibmparm_t,ProbClosures> ProbIB;

// ============================================================================
// update_geometry — apply rigid_trajectory to all geometry vertices
//
// NOTE: this deformable-geometry path is 2D-only (uses Polygon2D directly).
// In 3D builds it is elided so the TU still compiles; the rigid-body path
// in update_rigid_transforms is used instead.
// ============================================================================
#if (AMREX_SPACEDIM == 2)
inline void update_geometry(Real time, Vector<GeomType>& geom_a, int ngeom) {
  // Cache reference (t=0) vertex positions on first call
  static Vector<Polygon2D> ref_geom;
  if (ref_geom.empty()) {
    ref_geom.resize(ngeom);
    for (int i = 0; i < ngeom; i++) ref_geom[i] = geom_a[i];
    amrex::Print() << "[update_geometry] stored reference vertices ("
                   << ngeom << " geoms)\n";
  }

  RigidState s = rigid_trajectory(time);
  Real costh = std::cos(s.theta);
  Real sinth = std::sin(s.theta);
  Real cx = Motion::pitch_center_x;
  Real cy = Motion::pitch_center_y;

  // Transform: translate center, rotate about center, translate back + displacement
  for (int ig = 0; ig < ngeom; ig++) {
    int nv = static_cast<int>(ref_geom[ig].size());
    for (int iv = 0; iv < nv; iv++) {
      Real x0 = ref_geom[ig].vertex(iv)[0] - cx;
      Real y0 = ref_geom[ig].vertex(iv)[1] - cy;
      geom_a[ig].verts[iv][0] = cx + costh * x0 - sinth * y0 + s.dx;
      geom_a[ig].verts[iv][1] = cy + sinth * x0 + costh * y0 + s.dy;
    }
  }
}
#else
// 3D stub: deformable-geometry path is not implemented for TriMesh yet.
inline void update_geometry(Real /*time*/, Vector<GeomType>& /*geom_a*/, int /*ngeom*/) {
  amrex::Abort("update_geometry: 3D deformable path not implemented; use rigid transforms.");
}
#endif

// ============================================================================
// update_rigid_transforms — lightweight rigid-body path (no BVH rebuild)
//
// Computes rotation + translation and stores as a RigidTransform.
//
// IMPORTANT: we subtract the t=0 reference state so that the transform
// at t=0 is identity. This matches the geometry's actual initial position
// (loaded from file) and prevents a spurious "teleport" on the first step
// when rigid_trajectory(0) is nonzero (e.g. phase_y = π/2 Lissajous).
// ============================================================================
inline void update_rigid_transforms(Real time,
    Vector<RigidTransform>& transforms, int ngeom)
{
  // Relative state: subtract t=0 offset so the geometry starts at rest
  // in its loaded position, regardless of phase offsets.
  RigidState s0 = rigid_trajectory(0.0_rt);
  RigidState s  = rigid_trajectory(time);
  Real dx_rel    = s.dx    - s0.dx;
  Real dy_rel    = s.dy    - s0.dy;
  Real theta_rel = s.theta - s0.theta;

  Real costh = std::cos(theta_rel);
  Real sinth = std::sin(theta_rel);
  Real cx = Motion::pitch_center_x;
  Real cy = Motion::pitch_center_y;

  for (int i = 0; i < ngeom; i++) {
    RigidTransform& T = transforms[i];
    // 2D rotation matrix (relative to t=0 orientation)
    T.R[0][0] =  costh;  T.R[0][1] = -sinth;
    T.R[1][0] =  sinth;  T.R[1][1] =  costh;
    // Translation: (I - R) * pitch_center + relative displacement
    T.d[0] = cx - costh*cx + sinth*cy + dx_rel;
    T.d[1] = cy - sinth*cx - costh*cy + dy_rel;
  }
}

void inline inputs() {
  amrex::ParmParse ppm("motion");
  ppm.query("amp_x", Motion::amp_x);
  ppm.query("amp_y", Motion::amp_y);
  ppm.query("freq_x", Motion::freq_x);
  ppm.query("freq_y", Motion::freq_y);
  ppm.query("phase_y", Motion::phase_y);
  ppm.query("amp_theta", Motion::amp_theta);
  ppm.query("freq_theta", Motion::freq_theta);
  ppm.query("pitch_center_x", Motion::pitch_center_x);
  ppm.query("pitch_center_y", Motion::pitch_center_y);

  amrex::ParmParse ppib("ib");
  bool tmp_move = false;
  ppib.query("move", tmp_move);
  Motion::true_motion = tmp_move;

  amrex::Print() << " ****** Diamond Wedge — 3-DOF Oscillatory Motion ******* \n";
  amrex::Print() << " Translation: Ax=" << Motion::amp_x << " m, fx=" << Motion::freq_x
                 << " Hz, Ay=" << Motion::amp_y << " m, fy=" << Motion::freq_y
                 << " Hz, phase_y=" << Motion::phase_y << " rad\n";
  amrex::Print() << " Rotation: A_theta=" << Motion::amp_theta << " deg, f_theta="
                 << Motion::freq_theta << " Hz\n";
  amrex::Print() << " true_motion=" << Motion::true_motion << "\n";
  amrex::Print() << " ****************************************************** \n";
}

//////////////////////////// Initial conditions ////////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void prob_initdata (int i, int j, int k, amrex::Array4<amrex::Real> const& state,
      amrex::GeometryData const& geomdata, ProbClosures const& cls, ProbParm const& pp) {

  const Real* prob_lo = geomdata.ProbLo();
  const Real* dx      = geomdata.CellSize();
  Real x = prob_lo[0] + (i+0.5_rt)*dx[0];

  Real rho = pp.rho_oo;
  Real u0  = pp.u_oo;  // = 0 for quiescent atmosphere

  state(i, j, k, cls.URHO) = rho;
  state(i, j, k, cls.UMX)  = rho * u0;
  state(i, j, k, cls.UMY)  = 0.0;
  state(i, j, k, cls.UMZ)  = 0.0;
  state(i, j, k, cls.UET)  = rho * pp.eint_oo + 0.5_rt * rho * u0 * u0;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void user_tagging(int i, int j, int k, int nt, auto& tagfab, const auto &sdatafab,
                  const Array4<unsigned char>&ibfab, const auto& geomdata,
                  const ProbParm& pp, int level) {

  bool refine = false;

  if (level < 3) {
     bool foundGP = false;
     for (int ii = -1; ii <= 1; ii++) {
       for (int jj = -1; jj <= 1; jj++) {
#if (AMREX_SPACEDIM == 3)
         for (int kk = -1; kk <= 1; kk++) {
           foundGP = ibfab(i+ii,j+jj,k+kk,1) || foundGP;
         }
#else
         foundGP = ibfab(i+ii,j+jj,k,1) || foundGP;
#endif
       }
     }
     refine = foundGP || refine;
  }

  tagfab(i,j,k) = refine;
}


//////////////////////////// Boundary conditions ///////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const amrex::Real x[AMREX_SPACEDIM], amrex::Real dratio, const amrex::Real s_int[ProbClosures::NCONS],
         const amrex::Real s_refl[ProbClosures::NCONS], amrex::Real s_ext[ProbClosures::NCONS],
         const int idir, const int sgn, const amrex::Real time,
         amrex::GeometryData const& /*geomdata*/,  ProbClosures const& closures, ProbParm const& pp)
{
  const int URHO = ProbClosures::URHO;
  const int UMX  = ProbClosures::UMX;
  const int UMY  = ProbClosures::UMY;
  const int UET  = ProbClosures::UET;
  const int face = (idir+1)*sgn;

  Motion::sim_time = time;

  // Supersonic freestream Ma=2
  switch(face)
  {
    case  1:  // x-low: supersonic inflow
      s_ext[URHO] = pp.rho_oo;
      s_ext[UMX]  = pp.rho_oo * pp.u_oo;
      s_ext[UMY]  = 0.0_rt;
#if (AMREX_SPACEDIM == 3)
      s_ext[ProbClosures::UMZ] = 0.0;
#endif
      s_ext[UET]  = pp.rho_oo * pp.eint_oo + pp.kin_oo;
      break;
    case -1:  // x-high: supersonic outflow
      break;
    default:  // y-boundaries: outflow
      break;
  }
}

///////////////////////// IBM USER CLASS — moving wall BC
template < typename param, typename cls_t>
class ibm_user_t
{
  public:

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static void compute_surfIB(const Array1D<Real,0,AMREX_SPACEDIM-1>& xyz,
      const Array1D<Real,0,AMREX_SPACEDIM-1>& norm,
      const Array1D<Real,0,AMREX_SPACEDIM-1>& t1,
      const Array1D<Real,0,AMREX_SPACEDIM-1>& /*t2*/,
      Array2D<Real,0,ibm_eorder+1,0,cls_t::NPRIM-1>& q,
      int /*type_solid_bc*/, const cls_t* /*cls*/)
    {
      // Moving-wall no-slip BC from prescribed rigid-body kinematics.
      // Compute instantaneous velocity at surface point from rigid_trajectory.
      // V_wall = V_center + omega × (r - r_center)
      using namespace Motion;
      constexpr Real twopi = 2.0_rt * PI;
      Real t = sim_time;

      // Translational velocity
      Real uwx = amp_x * twopi * freq_x * std::cos(twopi * freq_x * t);
      Real uwy = amp_y * twopi * freq_y * std::cos(twopi * freq_y * t + phase_y);

      // Rotational contribution
      Real omega_deg = amp_theta * twopi * freq_theta * std::cos(twopi * freq_theta * t);
      Real om = deg2rad(omega_deg);
      Real rx = xyz(0) - pitch_center_x;
      Real ry = xyz(1) - pitch_center_y;
      uwx += -om * ry;
      uwy +=  om * rx;

      // Project wall velocity into local frame (n, t1)
      Real un  = uwx * norm(0) + uwy * norm(1);
      Real ut1 = uwx * t1(0) + uwy * t1(1);

      q(1,cls_t::QU) = un;
      q(1,cls_t::QV) = ut1;
#if (AMREX_SPACEDIM == 3)
      q(1,cls_t::QW) = 0.0_rt;
#endif

      q(1,cls_t::QPRES) = q(2,cls_t::QPRES);  // zero-gradient pressure
      q(1,cls_t::QT)    = param::Twall;         // isothermal wall
    }
};

}
#endif
