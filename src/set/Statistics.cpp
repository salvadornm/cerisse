// This file contains on-the-fly statistics functions.

#include <AMReX_ParmParse.H>

#include "CNS.h"

using namespace amrex;

//----------------------------------------------------------------------------------
// function to initialise statistic arrays to 0, is called per level
void CNS::setupStats() {

  // amrex::Print( ) << " Initialise STATS " << std::endl; 
  // amrex::Print( ) << " NSTATS = " << PROB::ProbClosures::NSTAT << std::endl; 
  // amrex:: Print( ) << " level = " << level << std::endl;
    
  MultiFab &S_stats = get_new_data(Stats_Type);
  auto const &data_stats = S_stats.arrays();
  amrex::ParallelFor(
    S_stats, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {
      for (int n=0; n < PROB::ProbClosures::NSTAT; n++)  {
        data_stats[box_no](i, j, k, n) =  0; 
      }                 
    });
}
//----------------------------------------------------------------------------------
// functon to compute and store mean values <U> and mean square <UU>, called per level
void CNS::computeStats(){
  // get arrays
  MultiFab &S_stats = get_new_data(Stats_Type);
  auto const &data_stats = S_stats.arrays();
  MultiFab &S = get_new_data(State_Type);
  auto const &data = S.arrays();

  Real dt = parent->dtLevel(level);
  
  // time_stats += dt; // obsolete
  Real time_aux = time_stat_level[level];
  
  // amrex::Print( ) << " Compute  STATS " << std::endl; 
  // amrex:: Print( ) << " level = " << level << std::endl;
  // amrex::Print( ) << " dt= " << dt <<  std::endl; 
  // amrex::Print( ) << " time_stat= " << time_stats <<  std::endl; 
  // amrex::Print( ) << " time_aux= " << time_aux <<  std::endl; 
  // amrex::Print( ) << " time_old= "  << time_stats-dt <<  std::endl; 
  // amrex::Print( ) << " counter_stat= " << counter_stats <<  std::endl; 
  

  // compute primitives ??? (need mfi)
  const PROB::ProbClosures& cls_h = *CNS::h_prob_closures;
  Array4<Real> const& state = statemf.array(mfi);
  const Box& bxg = mfi.growntilebox(cls_h.NGHOST);
  FArrayBox primf(bxg, cls_h.NPRIM, The_Async_Arena());
  Array4<Real> const& prims= primf.array();
  cls_h.cons2prims(mfi, state, prims); 
  ////

                        
  // storing stats
  amrex::ParallelFor(
    S_stats, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {       
    
    amrex::Real rho = data[box_no](i,j,k,PROB::ProbClosures::URHO);
    amrex::Real o_rho = 1.0/rho;

    // reset array  (so is sum (u*dt)
    for (int n=0; n< PROB::ProbClosures::NSTAT; n++){
      data_stats[box_no](i, j, k, n) *= (time_aux - dt);            
    }                

    if (PROB::ProbClosures::record_velocity > 0) { 
      amrex::Real vel[3]={0,0,0};         
      vel[0] = data[box_no](i,j,k,PROB::ProbClosures::UMX)*o_rho;
      vel[1] = data[box_no](i,j,k,PROB::ProbClosures::UMY)*o_rho;
      vel[2] = data[box_no](i,j,k,PROB::ProbClosures::UMZ)*o_rho;
       
      // U,V,W 
      for (int n=0; n< AMREX_SPACEDIM; n++){
        data_stats[box_no](i, j, k, n) += vel[n]*dt;            
      }                
      // UU,VV,WW        
      for (int n=0; n< AMREX_SPACEDIM ; n++) {
        data_stats[box_no](i, j, k, AMREX_SPACEDIM+n) += vel[n]*vel[n]*dt;            
      }       
      // UV,UW,VW        
#if (AMREX_SPACEDIM >= 2)          
      data_stats[box_no](i, j, k, 2*AMREX_SPACEDIM) += vel[0]*vel[1]*dt;
#endif
#if (AMREX_SPACEDIM == 3)
      data_stats[box_no](i, j, k, 2*AMREX_SPACEDIM+1) += vel[0]*vel[2]*dt;
      data_stats[box_no](i, j, k, 2*AMREX_SPACEDIM+2) += vel[1]*vel[2]*dt;
#endif
    }
    
    if (PROB::ProbClosures::record_PTrho > 0 ){

      //            
      Real P = prims(i, j, k, PROB::ProbClosures::QPRES);
      Real T = prims(i, j, k, PROB::ProbClosures::QT);
            
      // store P,T,rho
      data_stats[box_no](i, j, k, INDEX_THERM) += P*dt;
      data_stats[box_no](i, j, k, INDEX_THERM+1) += T*dt;
      data_stats[box_no](i, j, k, INDEX_THERM+2) += rho*dt;
      data_stats[box_no](i, j, k, INDEX_THERM+3) += P*P*dt;
      data_stats[box_no](i, j, k, INDEX_THERM+4) += T*T*dt;
      data_stats[box_no](i, j, k, INDEX_THERM+5) += rho*rho*dt;    
    }
    
    // store mean values  
    for (int n=0; n< PROB::ProbClosures::NSTAT; n++){
      data_stats[box_no](i, j, k, n) /= time_aux;            
    }           
    
  });
  
  // debugStats(10,96, 0);

}

//----------------------------------------------------------------------------------
// debugging function to access stat array
void CNS::debugStats(int is, int js, int ks) {
  
  MultiFab &S_stats = get_new_data(Stats_Type);
  auto const &data_stats = S_stats.arrays();

  printf(" *** DEBUG point >>>  i=%d j=%d k=%d \n",is,js,ks);
  amrex::ParallelFor(
    S_stats, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {
      
      if ((i==is) && (j==js))     {  
        printf(" box_no=%d \n ",box_no );
        for (int n=0; n< PROB::ProbClosures::NSTAT; n++) {      
        printf(" n=%d  DATAstat= %f \n ",n,data_stats[box_no](i, j, k, n) );
        }
      }         
    });    

  printf("  ***  \n ") ; 

}
