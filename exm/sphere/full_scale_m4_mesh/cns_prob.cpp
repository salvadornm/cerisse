// #include <CNS.H>
#include <AMReX_PROB_AMR_F.H>
#include <AMReX_ParmParse.H>
#include <cns_prob.H>

extern "C" {
    void amrex_probinit (const int* /*init*/,
                         const int* /*name*/,
                         const int* /*namelen*/,
                         const amrex_real* /*problo*/,
                         const amrex_real* /*probhi*/)
    {
        // could read parmparse parameters here

        // amrex::Gpu::htod_memcpy(CNS::d_prob_parm, CNS::h_prob_parm, sizeof(ProbParm));

        // tagging
    }
}


void cns_initdata (int i, int j, int k, amrex::Array4<amrex::Real> const& state,
              amrex::GeometryData const& geomdata, Parm const& parm, ProbParm const& prob_parm)
{
    using amrex::Real;

    const Real* prob_lo = geomdata.ProbLo();
    const Real* prob_hi = geomdata.ProbHi();
    const Real* dx      = geomdata.CellSize();

    Real x = prob_lo[0] + (i+Real(0.5))*dx[0];
    Real Pt, rhot, uxt;
    if (x < Real(0.5)*(prob_lo[0]+prob_hi[0])) {
        Pt   = prob_parm.p_l;
        rhot = prob_parm.rho_l;
        uxt  = prob_parm.u_l;
    } else {
        Pt   = prob_parm.p_r;
        rhot = prob_parm.rho_r;
        uxt  = prob_parm.u_r;
    }
    state(i,j,k,URHO ) = rhot;
    state(i,j,k,UMX  ) = rhot*uxt;
    state(i,j,k,UMY  ) = Real(0.0);
    state(i,j,k,UMZ  ) = Real(0.0);
    Real eint = Pt/(parm.eos_gamma-Real(1.0));
    state(i,j,k,UEINT) = eint;
    state(i,j,k,UEDEN) = eint + Real(0.5)*rhot*uxt*uxt;
    state(i,j,k,UTEMP) = Real(0.0);

    // amrex::Print() << i << " " << j << " " << k << " " << "\n";
    // amrex::Print() << prob_parm.p_r << "  "<< parm.eos_gamma << "\n";
    // exit(0);
}


void tagging(amrex::TagBoxArray& tags, amrex::MultiFab& sdata, int level, IBM::IBMultiFab* ibdata){

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  for (amrex::MFIter mfi(tags,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {
    const amrex::Box& bx = mfi.tilebox();
    auto const& tagfab = tags.array(mfi);
    auto const& sdatafab = sdata.array(mfi);

    if (ibdata!=nullptr) {
    auto const& ibdatafab = ibdata->array(mfi);
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      bool condition = false; 

      // tagging on ghs
      if (ibdatafab(i,j,k,1)) {
        condition = ibdatafab(i,j,k,1);
      }
      // or fluid neighbour
      else if (!ibdatafab(i,j,k,0)){
        for (int l=-1; l<=1 ; l+=2) {
        condition = condition || ibdatafab(i+l,j,k,0);
        condition = condition || ibdatafab(i,j+l,k,0);
        condition = condition || ibdatafab(i,j,k+l,0);
        }
      }
      // only tag front half of sphere on higher level
      if (level==1) {
        Real x=(0.5_rt + i)*IBM::ib.cellSizes[level][0];
        Real y=(0.5_rt + j)*IBM::ib.cellSizes[level][1];
        Real z=(0.5_rt + k)*IBM::ib.cellSizes[level][2];
        condition = condition && x<0.75_rt;
      }

      // store tag
      tagfab(i,j,k) = condition   ;
    });
    }
        // else  {
        //   // Temporary tagging on density in x only
        //   amrex::ParallelFor(bx,
        //   [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        //   {
        //   // amrex::Real drhox = amrex::Math::abs(sdatafab(i+1,j,k,0) - sdatafab(i-1,j,k,0))/sdatafab(i,j,k,0);
        //   // tagfab(i,j,k) = drhox > 0.5f;
        //   });
        // };
        // };


            // amrex::Print() << "i,j,k         " << i << " " << j << " " << k << " "<< std::endl;
            // amrex::Print() << "tag(i,j,k)    " << int(tagfab(i,j,k)) << std::endl;
            // amrex::Print() << "data(i,j,k,0) " << datafab(i,j,k,0) << std::endl;
            // amrex::Print() << "drhox         " << drhox << std::endl;
            // amrex::Print() << "------------- " << std::endl;

    }
}
