#ifndef UTILITY_H_
#define UTILITY_H_

#include <AMReX_AmrLevel.H>
#include <AMReX_FluxRegister.H>
#include <CNSconstants.h>

#ifdef USE_PELEPHYSICS
#include <PMFData.H>
#include <PMF.H>
#endif

#include <fstream>
#include <vector>
#include <iostream>

using namespace amrex;

// 
struct VelocityField {
    int Nx, Ny, Nz;
    std::vector<float> ux, uy, uz;
  }; 


// interface between PelePhysics<-> Cerisse as well as other utilities

class Utility{

  public:
  
  Utility()
  {};

  ~Utility()
  {};

  void test(){
    amrex::Print() << " ***  Utility test ***" << std::endl;
  }

// PMF: Premixed Flame Initilialization
#ifdef USE_PELEPHYSICS

  pele::physics::PMF::PmfData pmfData;

  // initialize PMF data
  void initPMF() {
    pmfData.initialize();
  }
#endif
  VelocityField turbData;

  // initialize Turbulence file
  void initTurbulenceFile(const std::string& filename){

    amrex::Print() << " ***  Reading Turbulence File ***" << std::endl;


    // read from input file

    turbData =  load_velocity_field(filename);

     
    // store database (where?)
  }   

  // returns velcoity given i,j,k in the mesh
  void get_velocity(const int i,const int j,const int k, 
      amrex::Real& u, amrex::Real& v, amrex::Real& w)
  {
    int Nx = turbData.Nx; int Ny = turbData.Ny;  
    u = turbData.ux[i + Nx*(j + Ny*k)];
    v = turbData.uy[i + Nx*(j + Ny*k)];
    w = turbData.uz[i + Nx*(j + Ny*k)];
  }


  //----------
  private:

  VelocityField load_velocity_field(const std::string& filename) {
    VelocityField field;
    std::ifstream file(filename, std::ios::binary);
    if (!file) throw std::runtime_error("File not found");

    // Read header
    file.read(reinterpret_cast<char*>(&field.Nx), sizeof(int));
    file.read(reinterpret_cast<char*>(&field.Ny), sizeof(int));
    file.read(reinterpret_cast<char*>(&field.Nz), sizeof(int));
    int N = field.Nx * field.Ny * field.Nz;

    // Read velocity components
    field.ux.resize(N);
    field.uy.resize(N);
    field.uz.resize(N);

    file.read(reinterpret_cast<char*>(field.ux.data()), N * sizeof(float));
    file.read(reinterpret_cast<char*>(field.uy.data()), N * sizeof(float));
    file.read(reinterpret_cast<char*>(field.uz.data()), N * sizeof(float));

    return field;
  }


};

#endif  
 

