#include "CNS.H"
// #include <AMReX_EBMultiFabUtil.H>

using namespace amrex;

void
CNS::restart (Amr& papa, std::istream& is, bool bReadSpecial)
{
  AmrLevel::restart(papa,is,bReadSpecial);

  if (do_reflux && level > 0) {
    flux_reg.define(grids, papa.boxArray(level-1),
                    dmap, papa.DistributionMap(level-1),
                    geom, papa.Geom(level-1),
                    papa.refRatio(level-1), level, LEN_STATE);
  }

  buildMetrics();
}

void
CNS::checkPoint (const std::string& dir, std::ostream& os, VisMF::How how, bool dump_old)
{
  AmrLevel::checkPoint(dir, os, how, dump_old);
}

void
CNS::writePlotFile (const std::string& dir, std::ostream& os, VisMF::How how)
{
  BL_PROFILE("CNS::writePlotFile()");
  AmrLevel::writePlotFile(dir, os, how);
}

void 
CNS::setPlotVariables ()
{
  // This will add everything into the plotfile by default
  amrex::AmrLevel::setPlotVariables();

  // So we remove the unwanted variables here
  amrex::ParmParse pp("cns");

  // Cost
  bool plot_cost = true;
  pp.query("plot_cost", plot_cost);
  if (!plot_cost) {
    for (int i = 0; i < desc_lst[Cost_Type].nComp(); i++) {
      amrex::Amr::deleteStatePlotVar(desc_lst[Cost_Type].name(i));
    }
  }

  // rho_omega and/or heat release
  if (!do_react) {
    for (int i = 0; i < desc_lst[Reactions_Type].nComp(); i++) {
      amrex::Amr::deleteStatePlotVar(desc_lst[Reactions_Type].name(i));
    }
  }

  bool plot_rho_omega = true;
  pp.query("plot_rho_omega", plot_rho_omega);
  if (do_react && !plot_rho_omega) {
    for (int nf = 0; nf <= NUM_FIELD; ++nf) {
    for (int i = 0; i < NUM_SPECIES + 1; i++) {
      amrex::Amr::deleteStatePlotVar(desc_lst[Reactions_Type].name(nf*NREACT + i));
    }
    }
  }

  bool plot_hrr = true;
  pp.query("update_heat_release", plot_hrr);
  if (do_react && !plot_hrr) {
    for (int nf = 0; nf <= NUM_FIELD; ++nf) {
      amrex::Amr::deleteStatePlotVar(desc_lst[Reactions_Type].name(nf*NREACT + NUM_SPECIES));
    }
  }

  // Transport coefficients
  if (!do_visc) {
    amrex::Amr::deleteDerivePlotVar("transport_coef");
  }

  // rhoY / Y / X
  bool plot_rhoy = true;
  pp.query("plot_rhoy", plot_rhoy);
  if (!plot_rhoy) {    
    for (int nf = 0; nf <= NUM_FIELD; ++nf) {
    for (int i = 0; i < NUM_SPECIES; i++) {
      amrex::Amr::deleteStatePlotVar(desc_lst[State_Type].name(nf*NVAR + UFS + i));
    }
    }
  }

  bool plot_massfrac = false;
  pp.query("plot_massfrac", plot_massfrac);
  if (plot_massfrac) {
    amrex::Amr::addDerivePlotVar("massfrac");
  } else {
    amrex::Amr::deleteDerivePlotVar("massfrac");
  }

  bool plot_molefrac = false;
  pp.query("plot_molefrac", plot_molefrac);
  if (plot_molefrac) {
    amrex::Amr::addDerivePlotVar("molefrac");
  } else {
    amrex::Amr::deleteDerivePlotVar("molefrac");
  }

  // Plot stochastic fields
  bool plot_fields = true;
  pp.query("plot_fields", plot_fields);
  if (!plot_fields) {
    for (int nf = 0; nf < NUM_FIELD; ++nf) {
      amrex::Amr::deleteStatePlotVar("density_Field" + std::to_string(nf));
      amrex::Amr::deleteStatePlotVar("xmom_Field" + std::to_string(nf));
      amrex::Amr::deleteStatePlotVar("ymom_Field" + std::to_string(nf));
      amrex::Amr::deleteStatePlotVar("zmom_Field" + std::to_string(nf));
      amrex::Amr::deleteStatePlotVar("rho_E_Field" + std::to_string(nf));
      for (int i = 0; i < NUM_SPECIES; ++i) {  
        amrex::Amr::deleteStatePlotVar("rho_omega_" + spec_names[i] + "_Field" + std::to_string(nf));
      }
      amrex::Amr::deleteStatePlotVar("heatRelease_Field" + std::to_string(nf));
    }
  }
}