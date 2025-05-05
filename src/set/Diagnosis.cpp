// This file contains on-the-fly diagnostics functions.

#include <AMReX_ParmParse.H>

#include "CNS.h"
#include <cassert>
#include <vector>

using namespace amrex;

// ======================== Time probe functionallity ========================
int CNS::time_probe_lev = 0;
int CNS::time_probe_int = 1;
Vector<std::string> CNS::time_probe_names;
Vector<Box> CNS::time_probe_boxes;

// Convert RealBox to Box. If RealBox does not align with grid, return a smaller
// box such that Box always > RealBox
Box realbox_to_box(std::vector<Real> const &rbox_lo,
                   std::vector<Real> const &rbox_hi, Geometry const &geom) {
  const Real *dx = geom.CellSize();
  const Real *prob_lo = geom.ProbLo();

  Box bx;
  for (int dir = 0; dir < amrex::SpaceDim; ++dir) {
    bx.setSmall(dir, ceil((rbox_lo[dir] - prob_lo[dir]) / dx[dir]));
    bx.setBig(dir, round((rbox_hi[dir] - prob_lo[dir]) / dx[dir]) - 1);
  }
  return bx;
}

void CNS::setupTimeProbe() {
//  assert(parent->NumDataLogs() > 0 &&
//         "No data logs file available for time probes!");

  ParmParse pp("cns");
  Vector<std::string> time_probes;
  pp.queryarr("time_probes", time_probes);
  const int num_probes = time_probes.size();
  if (num_probes == 0) return;

  pp.query("time_probe_lev", time_probe_lev);  // default 0
  pp.query("time_probe_int", time_probe_int);  // default 1

  for (int cnt = 0; cnt < num_probes; ++cnt) {
    ParmParse ppr(time_probes[cnt]);

    std::string field_name;
    ppr.get("field_name", field_name);
    time_probe_names.push_back(field_name);
    
    const Real *prob_lo = geom.ProbLo();
    const Real *prob_hi = geom.ProbHi();
    std::vector<Real> box_lo = {AMREX_D_DECL(prob_lo[0], prob_lo[1], prob_lo[2])};
    std::vector<Real> box_hi = {AMREX_D_DECL(prob_hi[0], prob_hi[1], prob_hi[2])};
    ppr.queryarr("box_lo", box_lo, 0, amrex::SpaceDim);
    ppr.queryarr("box_hi", box_hi, 0, amrex::SpaceDim);
    time_probe_boxes.push_back(realbox_to_box(box_lo, box_hi, geom));

    if (!time_probe_boxes[cnt].ok()) {      
      amrex::Abort("Invalid time probe box for " + field_name);
    }
  }

  // write header to file
  if (ParallelDescriptor::IOProcessor()) {
    const int log_index = 0;  // TODO: make this configurable?
    std::ostream &data_log = parent->DataLog(log_index);
    data_log << "time";
    for (int cnt = 0; cnt < num_probes; ++cnt) {
      data_log << ", " << time_probe_names[cnt] << "("
               << time_probe_boxes[cnt].smallEnd()
               << time_probe_boxes[cnt].bigEnd() << ")";
    }
    data_log << std::endl;
  }
}

// This needs:
//  int    time_probe_lev
//  int    verbose
//  int    time_probe_int
//  string time_probe_file
//  Vector<string> time_probe_names
//  Vector<Box>    time_probe_boxes
// It is user's responsibility to ensure ALL in_box are in time_probe_lev (by tagging for example),
//  as well as all in_box do NOT contain IB/EB covered cells.
void CNS::recordTimeProbe() {
  if (level != time_probe_lev) return;
  if (this->nStep() % time_probe_int != 0) return;

  if (verbose) {
    amrex::Print() << "... Processing time statistics\n";
  }

  const Real curtime = state[0].curTime();
  const int num_probes = time_probe_names.size(); // uncoment if C99 stantdard used
  //constexpr int num_probes = 1;
  //Real probe[num_probes] = {0.0};  // C99 standard (not compatible w/o tinkering flags)
  Real* probe = new Real[num_probes](); 
  MultiFab S(grids, dmap, h_prob_closures->NCONS, 1, MFInfo(), Factory());
  FillPatch(*this, S, 1, curtime, State_Type, 0, h_prob_closures->NCONS);

  for (int cnt = 0; cnt < num_probes; ++cnt) {
    const std::string &name = time_probe_names[cnt];
    const Box &in_box = time_probe_boxes[cnt];

    // if name is a state variable, just copy from state
    // if not, derive the box
    int index, scomp;
    if (isStateVariable(name, index, scomp)) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(S, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box bx = mfi.tilebox() & in_box;
        if (bx.ok()) {
          const auto &sarr = S[mfi].array();

          ReduceOps<ReduceOpSum> reduce_op;
          ReduceData<Real> reduce_data(reduce_op);
          using ReduceTuple = typename decltype(reduce_data)::Type;
          reduce_op.eval(
              bx, reduce_data,
              [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
                return sarr(i, j, k, scomp);
              });
          probe[cnt] += amrex::get<0>(reduce_data.value());
        }
      }
    } else if (const DeriveRec *rec = derive_lst.get(name)) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(S, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box bx = mfi.tilebox() & in_box;
        if (bx.ok()) {
          FArrayBox dfab(bx, 1);
          FArrayBox const &sfab = S[mfi];
          // we skip checking rec->derFuncFab() != nullptr
          rec->derFuncFab()(bx, dfab, 0, 1, sfab, geom, curtime, rec->getBC(),
                            level);
          const auto &darr = dfab.array();

          ReduceOps<ReduceOpSum> reduce_op;
          ReduceData<Real> reduce_data(reduce_op);
          using ReduceTuple = typename decltype(reduce_data)::Type;
          reduce_op.eval(
              bx, reduce_data,
              [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
                return darr(i, j, k, 0);
              });
          probe[cnt] += amrex::get<0>(reduce_data.value());
        }
      }
    } else {
      amrex::Abort("Unknown variable name in time_probe: " + name);
    }
  }

  // reduce over mpi
  ParallelDescriptor::ReduceRealSum(probe, num_probes,
                                    ParallelDescriptor::IOProcessorNumber());

  // write to file
  if (ParallelDescriptor::IOProcessor()) {
    const int log_index = 0;
    std::ostream &data_log = parent->DataLog(log_index);
    const int datprecision = 6;
    data_log << std::setprecision(datprecision) << curtime;
    for (int cnt = 0; cnt < num_probes; ++cnt) {
      data_log << ", " << std::setprecision(datprecision)
               << probe[cnt] / time_probe_boxes[cnt].numPts();
    }
    data_log << std::endl;
  }
}

// Note: output file name is specified by `amr.data_log`, the output should look like this
//  time, name1((0,0,0)(127,127,127)), name2((0,0,0)(127,127,127)), name3((64,64,64)(64,64,64)), ...
//  6.35576e-05, 0.0322707, 0.0965934, 0.000581789, 18.5371
//  0.000127115, 0.0322706, 0.0965934, 0.000581789, 18.5371
//  0.000190673, 0.0322706, 0.0965934, 0.000581789, 18.5371
//  0.00025423, 0.0322706, 0.0965935, 0.000581789, 18.5371
//  0.000317786, 0.0322706, 0.0965936, 0.000581789, 18.5371

// ===========================================================================