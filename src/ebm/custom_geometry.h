#ifndef CNS_CUSTOM_GEOMETRY_H_
#define CNS_CUSTOM_GEOMETRY_H_

#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <Factory2.H>

using namespace amrex;

/**
 * @brief Base class to register custom geometries.
 */
class CustomGeometry : public cerisse::physics::Factory<CustomGeometry>
{
public:
  static const std::string base_identifier() { return "CustomGeometry"; }

  ~CustomGeometry() override = default;

  virtual void build(const Geometry& geom, const int max_coarsening_level) = 0;
};

//
// Declare custom geometries
// See PelePhysics Factory.H for how to define a subclass

#if AMREX_SPACEDIM == 3
class Combustor : public CustomGeometry::Register<Combustor>
{
public:
  static const std::string identifier() { return "combustor"; }

  void build(const Geometry& geom, const int max_coarsening_level) override;
};

class ConvergingNozzle : public CustomGeometry::Register<ConvergingNozzle>
{
public:
  static std::string identifier() { return "converging-nozzle"; }

  void build(const Geometry& geom, const int max_coarsening_level) override;
};
#endif

class Triangles : public CustomGeometry::Register<Triangles>
{
public:
  static std::string identifier() { return "triangles"; }

  void build(const Geometry& geom, const int max_coarsening_level) override;
};

//
// Define methods to build custom geometries in custom_geometry.cpp

class Custom : public CustomGeometry::Register<Custom>
{
public:
  static const std::string identifier() { return "custom"; }

  void build(const Geometry& geom, const int max_coarsening_level) override;
};


#endif