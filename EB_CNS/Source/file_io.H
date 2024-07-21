#ifndef CNS_FILE_IO_H_
#define CNS_FILE_IO_H_

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <fstream>
#include <sstream>

AMREX_FORCE_INLINE
std::string read_file(std::ifstream& in)
{
  return static_cast<std::stringstream const&>(std::stringstream() << in.rdbuf())
    .str();
}

void read_binary(const std::string& iname, const size_t nx, const size_t ny,
                 const size_t nz, const size_t ncol,
                 amrex::Vector<double>& data /*needs to be double*/);

void read_csv(const std::string& iname, const size_t nx, const size_t ny,
              const size_t nz, amrex::Vector<amrex::Real>& data);

#endif