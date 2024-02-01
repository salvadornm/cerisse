#ifndef THERMODYNAMICS_H_
#define THERMODYNAMICS_H_

using namespace amrex;

template <typename idx_t>
class calorifically_perfect_gas_t {
 protected:
 idx_t idx;

 public:
  Real gamma = 1.40;   // ratio of specific heats
  Real mw = 28.96e-3;  // mean molecular weight air kg/mol

  Real Ru = Real(8.31451);
  Real cv = Ru / (mw * (gamma - Real(1.0)));
  Real cp = gamma * Ru / (mw * (gamma - Real(1.0)));
  Real Rspec = Ru / mw;

  Real sos(Real& T){return std::sqrt(this->gamma * this->Rspec * T);}

  // void prims2cons(i,j,k,){};

  void prims2char(){};


  // void prims2flux(int& i, int& j, int& k, const GpuArray<Real,NPRIM>& prims, GpuArray<Real,NCONS>& fluxes, const GpuArray<int, 3>& vdir) {
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  void prims2fluxes(int& i, int& j, int& k, const Array4<Real>& prims, Array4<Real>& fluxes, const GpuArray<int, 3>& vdir) {

    Real rho = prims(i,j,k,idx.QRHO);
    Real ux  = prims(i,j,k,idx.QU);
    Real uy  = prims(i,j,k,idx.QV);
    Real uz  = prims(i,j,k,idx.QW);
    Real P   = prims(i,j,k,idx.QPRES);
    Real udir = ux*vdir[0] + uy*vdir[1] + uz*vdir[2];

    Real ekin  = Real(0.5) *(ux*ux + uy*uy + uz*uz);
    Real rhoet = rho*(this->cp*prims(i,j,k,idx.QT) + ekin) ;

    fluxes(i,j,k,idx.URHO) = rho*udir;
    fluxes(i,j,k,idx.UMX)  = rho*ux*udir + P*vdir[0];
    fluxes(i,j,k,idx.UMY)  = rho*uy*udir + P*vdir[1];
    fluxes(i,j,k,idx.UMZ)  = rho*uz*udir + P*vdir[2];
    fluxes(i,j,k,idx.UET)  = (rhoet + P) * udir;
  };

  // can move this to closures derived tyoe (closures_dt)
  // prims to cons
  // - We want to call it from thermodynamics class
  // - cls is stored on cpu and gpu
  // TODO: remove ParallelFor from here. Keep closures local
  void inline cons2prims(const MFIter& mfi, const Array4<Real>& cons,
                         const Array4<Real>& prims) const {
    const Box& bxg = mfi.growntilebox(NGHOST);

    amrex::ParallelFor(bxg, [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) {
      Real rho = cons(i, j, k, idx.URHO);
      // Print() << "cons2prim"<< i << j << k << rho << std::endl;
      rho = max(1e-40, rho);
      Real rhoinv = Real(1.0) / rho;
      Real ux = cons(i, j, k, idx.UMX) * rhoinv;
      Real uy = cons(i, j, k, idx.UMY) * rhoinv;
      Real uz = cons(i, j, k, idx.UMZ) * rhoinv;
      Real rhoke = Real(0.5) * rho * (ux * ux + uy * uy + uz * uz);
      Real rhoei = (cons(i, j, k, idx.UET) - rhoke);
      Real p = (this->gamma - Real(1.0)) * rhoei;

      prims(i, j, k, idx.QRHO) = rho;
      prims(i, j, k, idx.QU) = ux;
      prims(i, j, k, idx.QV) = uy;
      prims(i, j, k, idx.QW) = uz;
      prims(i, j, k, idx.QPRES) = p;
      prims(i, j, k, idx.QT) = p / (rho * (this->Rspec));
    });
  }
};

class calorifically_perfect_gas_nasg_liquid_t
{
private:
  /* data */
  // gamma_a =;
  // gamma_l =;

public:
// state

Real inline energy (Real p, Real rho) {
  Real eint=0;
return eint;
}

// pressure (eint,rho)

};




#endif
