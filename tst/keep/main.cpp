#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BoxIterator.H>
#include <cmath>
#include <limits>
#include <iostream>

#ifndef KEEP_TEST_SPECIES
#define KEEP_TEST_SPECIES 3
#endif
static constexpr int NUM_SPECIES = KEEP_TEST_SPECIES;
using namespace amrex;
#include <Keep.h>
#include <CentralKEEP_SPEC.h>

// Composition-dependent EOS, deliberately without a cv member. This catches
// accidentally reverting to the calorically-perfect-gas-only KEEP flux.
struct Closure {
  static constexpr int QRHO=0, QU=1, QV=2, QW=3, QT=4, QPRES=5, QFS=6;
  static constexpr int QC=QFS+NUM_SPECIES, NPRIM=QC+1;
  static constexpr int UMX=0, UMY=1, UMZ=2, UET=3, UFS=4, URHO=UFS;
  static constexpr int NCONS=UFS+NUM_SPECIES, NGHOST=3;
  void TYR2E(Real T, const Real* Y, Real, Real& e) const {
    e=0;
    for(int n=0;n<NUM_SPECIES;++n) e += Real(n+2)*Y[n]*T;
  }
};
template<int N, bool D=false> struct Params {
  static constexpr int order=N;
  static constexpr bool dissipation=D;
  static constexpr Real C2skew=.1, C4skew=.016;
};
void check(bool ok, const char* message) {
  if(!ok) amrex::Abort(message);
}
void close(Real a, Real b, const char* message) {
  check(std::isfinite(a) && std::isfinite(b) &&
        std::abs(a-b)<Real(2e-12)*amrex::max(Real(1),std::abs(b)),message);
}

void fill(const Box& bx, Array4<Real> q, Array4<Real> u, int dir, bool uniform) {
  for (BoxIterator it(bx); it.ok(); ++it) {
    const IntVect iv=it();
    const Real s=uniform ? 0 : Real(iv[dir]);
    const Real rho=1.2+.02*s;
    Real Y[NUM_SPECIES], total=0;
    for(int n=0;n<NUM_SPECIES;++n) {Y[n]=Real(n+1)+.01*s*(n+2); total+=Y[n];}
    for(int n=0;n<NUM_SPECIES;++n) {Y[n]/=total; q(iv,Closure::QFS+n)=Y[n];}
    const Real T=2+.03*s;
    Real e; Closure{}.TYR2E(T,Y,rho,e);
    q(iv,Closure::QRHO)=rho;
    q(iv,Closure::QT)=T;
    q(iv,Closure::QPRES)=.4*rho*e + .001*s*s;
    q(iv,Closure::QC)=2;
    Real ke=0;
    for(int d=0;d<3;++d) {
      const Real v=.7+.02*s;
      q(iv,Closure::QU+d)=v;
      u(iv,Closure::UMX+d)=rho*v;
      ke+=.5*v*v;
    }
    u(iv,Closure::UET)=rho*(e+ke);
    for(int n=0;n<NUM_SPECIES;++n) u(iv,Closure::UFS+n)=rho*Y[n];
  }
}

template<int N> void test_flux() {
  const Box bx(IntVect(-4),IntVect(4));
  FArrayBox qfab(bx,Closure::NPRIM), ufab(bx,Closure::NCONS);
  FArrayBox fnew(bx,Closure::NCONS), fold(bx,Closure::NCONS), fwall(bx,Closure::NCONS);
  BaseFab<uint8_t> mask(bx,2);
  auto q=qfab.array(),u=ufab.array(),f=fnew.array(),g=fold.array(),h=fwall.array();
  auto mark=mask.array();
  Closure cls;
  keep_t<Params<N>,Closure> keep;
  keep_t<Params<2>,Closure> two;
  keep_t<Params<N,true>,Closure> damped;
  keep_euler_spec_t<false,false,N,Closure> old;
  Array1D<Real,0,2> coeff = N==2 ? Array1D<Real,0,2>{1,0,0} :
      (N==4 ? Array1D<Real,0,2>{4./3.,-1./6.,0} : Array1D<Real,0,2>{1.5,-.3,1./30.});
  Real wall_mass[AMREX_SPACEDIM];
  for(int dir=0;dir<AMREX_SPACEDIM;++dir) {
    IntVect left(0); left[dir]=-1;
    fill(bx,q,u,dir,false); mask.setVal(0);
    keep.flux_dir(0,0,0,dir,q,u,f,&cls);
    old.flux_dir(0,0,0,dir,coeff,q,g,&cls);
    keep.flux_dir_ibm(0,0,0,dir,q,u,h,&cls,mark);
    for(int n=0;n<Closure::NCONS;++n) {
      close(f(0,0,0,n),g(0,0,0,n),"KEEP differs from legacy split flux");
      close(h(0,0,0,n),f(0,0,0,n),"regular marked flux differs");
    }
    // Independently sum the split pair mass flux to check species conservation.
    Real expected_mass=0, sum=0;
    for(int l=1;l<=N/2;++l) for(int m=0;m<l;++m) {
      IntVect a(0),b(0);a[dir]=m;b[dir]=m-l;
      expected_mass+=coeff(l-1)*.25*(q(a,Closure::QRHO)+q(b,Closure::QRHO))*
          (q(a,Closure::QU+dir)+q(b,Closure::QU+dir));
    }
    for(int n=0;n<NUM_SPECIES;++n) sum+=f(0,0,0,Closure::UFS+n);
    close(sum,expected_mass,"species fluxes do not sum to split mass flux");

    // Both boundary modes must keep the high-order dissipative interior flux.
    damped.flux_dir(0,0,0,dir,q,u,f,&cls);
    damped.flux_dir_ibm(0,0,0,dir,q,u,h,&cls,mark);
    for(int n=0;n<Closure::NCONS;++n) close(f(0,0,0,n),h(0,0,0,n),"regular dissipative marker mismatch");

    // EB cut cell and IBM ghost point: fallback is the SAME second-order KEEP.
    for(int solid : {0,1}) {
      mask.setVal(0); mark(left,0)=solid; mark(left,1)=1;
      keep.flux_dir_ibm(0,0,0,dir,q,u,h,&cls,mark);
      two.flux_dir(0,0,0,dir,q,u,g,&cls);
      for(int n=0;n<Closure::NCONS;++n) close(h(0,0,0,n),g(0,0,0,n),"wall fallback differs from two-point KEEP");
    }

    // Poison all non-adjacent stencil states as covered solid. Neither the
    // EOS flux nor dissipation may touch them, including for y and z faces.
    mask.setVal(0);
    for(int l=-3;l<=2;++l) if(l!=-1 && l!=0) {
      IntVect iv(0);iv[dir]=l;mark(iv,0)=1;
      for(int n=0;n<Closure::NPRIM;++n) q(iv,n)=std::numeric_limits<Real>::quiet_NaN();
      for(int n=0;n<Closure::NCONS;++n) u(iv,n)=std::numeric_limits<Real>::quiet_NaN();
    }
    keep.flux_dir_ibm(0,0,0,dir,q,u,h,&cls,mark);
    for(int n=0;n<Closure::NCONS;++n) close(h(0,0,0,n),g(0,0,0,n),"covered remote state entered KEEP stencil");
    damped.flux_dir_ibm(0,0,0,dir,q,u,h,&cls,mark);
    for(int n=0;n<Closure::NCONS;++n) check(std::isfinite(h(0,0,0,n)),"wall dissipation accessed covered state");
    wall_mass[dir]=h(0,0,0,Closure::UFS);

    // Interior solid face must clear old flux without evaluating poisoned EOS.
    mark(IntVect(0),0)=1;mark(left,0)=1;
    fwall.setVal(99);
    damped.flux_dir_ibm(0,0,0,dir,q,u,h,&cls,mark);
    for(int n=0;n<Closure::NCONS;++n) check(h(0,0,0,n)==0,"solid face not zeroed");

    // Constant-state consistency against the physical Euler flux, including e(Y,T).
    fill(bx,q,u,dir,true);
    damped.flux_dir(0,0,0,dir,q,u,f,&cls);
    const Real rho=q(0,0,0,Closure::QRHO),v=q(0,0,0,Closure::QU+dir),p=q(0,0,0,Closure::QPRES);
    for(int d=0;d<3;++d) close(f(0,0,0,Closure::UMX+d),rho*v*.7+(d==dir?p:0),"constant momentum flux");
    close(f(0,0,0,Closure::UET),(u(0,0,0,Closure::UET)+p)*v,"constant total-energy flux");
    for(int n=0;n<NUM_SPECIES;++n) close(f(0,0,0,Closure::UFS+n),v*u(0,0,0,Closure::UFS+n),"constant species flux");
  }
  for(int dir=1;dir<AMREX_SPACEDIM;++dir) close(wall_mass[dir],wall_mass[0],"wall sensor is direction-dependent");
}

// Exercise the real MFIter/face-box entry point selected by the EB/IBM macros.
template<bool D> void test_entry_point() {
  const Box valid(IntVect(0),IntVect(3));
  BoxArray ba(valid);DistributionMapping dm(ba);
  MultiFab q(ba,dm,Closure::NPRIM,3),u(ba,dm,Closure::NCONS,3);
  Closure cls;keep_t<Params<6,D>,Closure> method;
  for(MFIter mfi(q);mfi.isValid();++mfi) {
    fill(mfi.fabbox(),q.array(mfi),u.array(mfi),0,true);
    BaseFab<uint8_t> markers(mfi.fabbox(),2);markers.setVal(0);
    std::array<FArrayBox,AMREX_SPACEDIM> storage;
    std::array<FArrayBox*,AMREX_SPACEDIM> ptrs;
    for(int dir=0;dir<AMREX_SPACEDIM;++dir) {
      storage[dir].resize(surroundingNodes(valid,dir),Closure::NCONS);
      storage[dir].setVal(std::numeric_limits<Real>::quiet_NaN());
      ptrs[dir]=&storage[dir];
    }
#if (AMREX_USE_GPIBM || CNS_USE_EB)
    method.eflux_ibm(Geometry{},mfi,q.array(mfi),ptrs,u.array(mfi),&cls,markers.array());
#else
    method.eflux(Geometry{},mfi,q.array(mfi),ptrs,u.array(mfi),&cls);
#endif
    Gpu::streamSynchronize();
    for(int dir=0;dir<AMREX_SPACEDIM;++dir) {
      const auto flux=storage[dir].array();
      for(BoxIterator it(storage[dir].box());it.ok();++it) {
        Real mass=0;for(int n=0;n<NUM_SPECIES;++n)mass+=flux(it(),Closure::UFS+n);
        close(mass,Real(1.2*.7),"entry point did not fill every normal face");
      }
    }
  }
}

int main(int argc,char** argv) {
  amrex::Initialize(argc,argv);
  {
    test_flux<2>();test_flux<4>();test_flux<6>();
    test_entry_point<false>();test_entry_point<true>();
    amrex::Print()<<"PASS: KEEP orders 2/4/6; legacy equivalence, species conservation, EOS energy,\n"
                  <<"constant states, EB/IBM wall fallback, solid masks, directional dissipation, face loops.\n";
  }
  amrex::Finalize();
}
