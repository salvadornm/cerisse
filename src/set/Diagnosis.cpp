// This file contains on-the-fly diagnostics functions.

#include <AMReX_ParmParse.H>

#include "CNS.h"
#include <cassert>
#include <limits>
#include <vector>

using namespace amrex;

// ======================== Time probe functionallity ========================
int CNS::time_probe_lev = 0;
int CNS::time_probe_int = 1;

Vector<std::string> CNS::time_probe_names;
Vector<std::string> CNS::time_probe_reductions;
Vector<Box> CNS::time_probe_boxes;
Vector<int> CNS::time_probe_components;

// Convert RealBox to Box. If RealBox does not align with grid, return a smaller
// box such that Box always > RealBox
Box realbox_to_box(std::vector<Real> const &rbox_lo,
                   std::vector<Real> const &rbox_hi, Geometry const &geom) {
  const Real *dx = geom.CellSize();
  const Real *prob_lo = geom.ProbLo();

  Box bx;
  
  for (int dir = 0; dir < amrex::SpaceDim; ++dir) {
    //  bx.setSmall(dir, ceil((rbox_lo[dir] - prob_lo[dir]) / dx[dir]));
    //  bx.setBig(dir, round((rbox_hi[dir] - prob_lo[dir]) / dx[dir]) - 1);
    const Real ilo = (rbox_lo[dir] - prob_lo[dir]) / dx[dir] - Real(0.5);
    const Real ihi = (rbox_hi[dir] - prob_lo[dir]) / dx[dir] - Real(0.5);
    bx.setSmall(dir, static_cast<int>(std::ceil(ilo)));
    bx.setBig  (dir, static_cast<int>(std::floor(ihi)));
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

  if (time_probe_lev < 0 || time_probe_lev > parent->maxLevel()) {
    amrex::Abort("time_probe_lev must be between 0 and amr.max_level");
  }

  // Probes are recorded at time_probe_lev; use that level's Geometry so the
  // resulting Box indices live in the same index space as tiles at recording.
  const Geometry& probe_geom = parent->Geom(time_probe_lev);

  // Clear in case of restart (setupTimeProbe is called again in post_restart).
  time_probe_names.clear();
  time_probe_reductions.clear();
  time_probe_boxes.clear();
  time_probe_components.clear();


  for (int cnt = 0; cnt < num_probes; ++cnt) {
    ParmParse ppr(time_probes[cnt]);

    std::string field_name;
    ppr.get("field_name", field_name);
    time_probe_names.push_back(field_name);

    int component = 0;
    ppr.query("component", component);
    time_probe_components.push_back(component);

    // Optional: average/mean/avg (default), max/maximum, min/minimum
    std::string reduction = "average";
    ppr.query("reduction", reduction);
    if (reduction == "avg" || reduction == "mean") reduction = "average";
    if (reduction == "maximum") reduction = "max";
    if (reduction == "minimum") reduction = "min";
    if (reduction != "average" && reduction != "max" && reduction != "min") {
      amrex::Abort("Invalid time probe reduction for " + time_probes[cnt] +
                   ": use average, max, or min");
    }
    time_probe_reductions.push_back(reduction);

    const Real *prob_lo = probe_geom.ProbLo();
    const Real *prob_hi = probe_geom.ProbHi();
    std::vector<Real> box_lo = {AMREX_D_DECL(prob_lo[0], prob_lo[1], prob_lo[2])};
    std::vector<Real> box_hi = {AMREX_D_DECL(prob_hi[0], prob_hi[1], prob_hi[2])};
    ppr.queryarr("box_lo", box_lo, 0, amrex::SpaceDim);
    ppr.queryarr("box_hi", box_hi, 0, amrex::SpaceDim);
    time_probe_boxes.push_back(realbox_to_box(box_lo, box_hi, probe_geom));

    if (!time_probe_boxes[cnt].ok()) {
    amrex::Abort(
        "Invalid time probe box for probe '" + time_probes[cnt] +
        "', field '" + field_name + "'");
    }   

  }

  // write header to file
  if (ParallelDescriptor::IOProcessor()) {
    const int log_index = 0;  // TODO: make this configurable?
    std::ostream &data_log = parent->DataLog(log_index);
    data_log << "time";
    for (int cnt = 0; cnt < num_probes; ++cnt) {
      data_log << ", " << time_probe_names[cnt] << "_"
               << time_probe_reductions[cnt] << "("
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
  const int num_probes = time_probe_names.size();
  if (num_probes == 0) return;

  Vector<Real> probe(num_probes, Real(0.0));
  for (int cnt = 0; cnt < num_probes; ++cnt) {
    if (time_probe_reductions[cnt] == "max") {
      probe[cnt] = -std::numeric_limits<Real>::max();
    } else if (time_probe_reductions[cnt] == "min") {
      probe[cnt] =  std::numeric_limits<Real>::max();
    }
  }

  MultiFab S(grids, dmap, h_prob_closures->NCONS, 1, MFInfo(), Factory());
  FillPatch(*this, S, 1, curtime, State_Type, 0, h_prob_closures->NCONS);

  auto accumulate_probe = [&] (int cnt, Box const& bx, Array4<Real const> const& arr,
                               int comp) {
    const std::string& reduction = time_probe_reductions[cnt];

    if (reduction == "average") {
      ReduceOps<ReduceOpSum> reduce_op;
      ReduceData<Real> reduce_data(reduce_op);
      using ReduceTuple = typename decltype(reduce_data)::Type;
      reduce_op.eval(bx, reduce_data,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
          return arr(i, j, k, comp);
        });
      probe[cnt] += amrex::get<0>(reduce_data.value());
    } else if (reduction == "max") {
      ReduceOps<ReduceOpMax> reduce_op;
      ReduceData<Real> reduce_data(reduce_op);
      using ReduceTuple = typename decltype(reduce_data)::Type;
      reduce_op.eval(bx, reduce_data,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
          return arr(i, j, k, comp);
        });
      probe[cnt] = amrex::max(probe[cnt], amrex::get<0>(reduce_data.value()));
    } else { // min
      ReduceOps<ReduceOpMin> reduce_op;
      ReduceData<Real> reduce_data(reduce_op);
      using ReduceTuple = typename decltype(reduce_data)::Type;
      reduce_op.eval(bx, reduce_data,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
          return arr(i, j, k, comp);
        });
      probe[cnt] = amrex::min(probe[cnt], amrex::get<0>(reduce_data.value()));
    }
  };

  for (int cnt = 0; cnt < num_probes; ++cnt) {
    const std::string &name = time_probe_names[cnt];
    const Box &in_box = time_probe_boxes[cnt];

    int index, scomp;
    if (isStateVariable(name, index, scomp)) {
      for (MFIter mfi(S, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box bx = mfi.tilebox() & in_box;
        if (bx.ok()) {
          accumulate_probe(cnt, bx, S[mfi].const_array(), scomp);
        }
      }
    } 
    // old
    // else if (const DeriveRec *rec = derive_lst.get(name)) {
    //   for (MFIter mfi(S, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    //     const Box bx = mfi.tilebox() & in_box;
    //     if (bx.ok()) {
    //       FArrayBox dfab(bx, 1);
    //       FArrayBox const &sfab = S[mfi];
    //       rec->derFuncFab()(bx, dfab, 0, 1, sfab, geom, curtime, rec->getBC(),
    //                         level);
    //       accumulate_probe(cnt, bx, dfab.const_array(), 0);
    //     }
    //   }
    // } 
    // new
    else if (const DeriveRec *rec = derive_lst.get(name)) {
      const int nderive = rec->numDerive();
      const int selected_comp = time_probe_components[cnt];
      if (selected_comp < 0 || selected_comp >= nderive) {
      amrex::Abort(
          "Invalid component " + std::to_string(selected_comp) +
        " for derived variable " + name +
        ", which has " + std::to_string(nderive) + " components");
      }

      for (MFIter mfi(S, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box bx = mfi.tilebox() & in_box;

        if (bx.ok()) {
          /* The derive routine may write all components, so allocate all of them.*/
          FArrayBox dfab(bx, nderive);
          FArrayBox const& sfab = S[mfi];
          rec->derFuncFab()(bx, dfab, 0, nderive, sfab, geom, curtime, rec->getBC(), level);
          accumulate_probe(cnt,bx, dfab.const_array(), selected_comp);
        }
      }
    }  
    //
    else {
      amrex::Abort("Unknown variable name in time_probe: " + name);
    }
  }

  // Reduce over MPI. Average probes reduce sums; max/min probes reduce extrema.
  Vector<int> avg_ids, max_ids, min_ids;
  Vector<Real> avg_vals, max_vals, min_vals;
  for (int cnt = 0; cnt < num_probes; ++cnt) {
    if (time_probe_reductions[cnt] == "average") {
      avg_ids.push_back(cnt); avg_vals.push_back(probe[cnt]);
    } else if (time_probe_reductions[cnt] == "max") {
      max_ids.push_back(cnt); max_vals.push_back(probe[cnt]);
    } else {
      min_ids.push_back(cnt); min_vals.push_back(probe[cnt]);
    }
  }

  if (!avg_vals.empty()) {
    ParallelDescriptor::ReduceRealSum(avg_vals.data(), avg_vals.size(),
                                      ParallelDescriptor::IOProcessorNumber());
    for (int n = 0; n < avg_ids.size(); ++n) probe[avg_ids[n]] = avg_vals[n];
  }
  if (!max_vals.empty()) {
    ParallelDescriptor::ReduceRealMax(max_vals.data(), max_vals.size(),
                                      ParallelDescriptor::IOProcessorNumber());
    for (int n = 0; n < max_ids.size(); ++n) probe[max_ids[n]] = max_vals[n];
  }
  if (!min_vals.empty()) {
    ParallelDescriptor::ReduceRealMin(min_vals.data(), min_vals.size(),
                                      ParallelDescriptor::IOProcessorNumber());
    for (int n = 0; n < min_ids.size(); ++n) probe[min_ids[n]] = min_vals[n];
  }

  if (ParallelDescriptor::IOProcessor()) {
    const int log_index = 0;
    std::ostream &data_log = parent->DataLog(log_index);
    const int datprecision = 6;
    data_log << std::setprecision(datprecision) << curtime;
    for (int cnt = 0; cnt < num_probes; ++cnt) {
      Real value = probe[cnt];
      if (time_probe_reductions[cnt] == "average") {
        value /= time_probe_boxes[cnt].numPts();
      }
      data_log << ", " << std::setprecision(datprecision) << value;
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