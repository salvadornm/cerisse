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





void user_tagging(amrex::TagBoxArray& tags, amrex::MultiFab& sdata, int level, IBM::IBMultiFab* ibdata){

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
            // tagging on ghs
            tagfab(i,j,k) = ibdatafab(i,j,k,1);
            // or fluid neighbour
            bool condition = false; 
            for (int l=-1; l<=1 ; l+=2) {
              condition = condition || ibdatafab(i+l,j,k,0);
              condition = condition || ibdatafab(i,j+l,k,0);
              condition = condition || ibdatafab(i,j,k+l,0);
            }
            condition = condition && !ibdatafab(i,j,k,0);
            tagfab(i,j,k) = tagfab(i,j,k) || condition   ;
          });}

        else  {
          // Temporary tagging on density in x only
          amrex::ParallelFor(bx,
          [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
          // amrex::Real drhox = amrex::Math::abs(sdatafab(i+1,j,k,0) - sdatafab(i-1,j,k,0))/sdatafab(i,j,k,0);
          // tagfab(i,j,k) = drhox > 0.5f;
          });
        };
        };


            // amrex::Print() << "i,j,k         " << i << " " << j << " " << k << " "<< std::endl;
            // amrex::Print() << "tag(i,j,k)    " << int(tagfab(i,j,k)) << std::endl;
            // amrex::Print() << "data(i,j,k,0) " << datafab(i,j,k,0) << std::endl;
            // amrex::Print() << "drhox         " << drhox << std::endl;
            // amrex::Print() << "------------- " << std::endl;

    }
