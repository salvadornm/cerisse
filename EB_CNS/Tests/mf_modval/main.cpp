#include <AMReX.H>
#include <AMReX_Amr.H>
#include <AMReX_ParallelDescriptor.H>

#include "CNS.H"

amrex::LevelBld* getLevelBld();
// #if CNS_USE_EB
// void initialize_EB2(const Geometry& geom, const int required_level,
//                     const int max_level);
// #endif

class TwoAmr
{
  class HelperAmr : public Amr
  {
    friend TwoAmr;

  public:
    template <typename... T>
    HelperAmr(T&&... args) : Amr(std::forward<T>(args)...)
    {
    }
  };

public:
  TwoAmr()
  {
    const int coord = 0;
    ParmParse pp("manip");

    amrex::Print() << ">>> Init amr1 ...\n";
    auto [rb_in, max_level_in, n_cell_in] = get_amr1();
    HelperAmr amr1(&rb_in, max_level_in, n_cell_in, coord, getLevelBld());
    amrex::Print() << "amr1.Geom(): \n";
    for (int lev = 0; lev <= amr1.maxLevel(); ++lev) {
      amrex::Print() << lev << ": " << amr1.Geom(lev).ProbDomain()
                     << amr1.Geom(lev).Domain() << amr1.Geom(lev).CellSize(0)
                     << "\n";
    }
    amrex::Print() << ">>> Reading chkfile " << amr1.restart_chkfile
                   << " to amr1 ...\n";
    amr1.restart(amr1.restart_chkfile);

    amrex::Print() << ">>> Init amr2 ...\n";
    tie(rb_in, max_level_in, n_cell_in) = get_amr2();
    int fill_level = 0; pp.query("fill_level", fill_level);
    HelperAmr amr2(&rb_in, fill_level, n_cell_in, coord, getLevelBld());
    amrex::Print() << "amr2.Geom(): \n";
    for (int lev = 0; lev <= amr2.maxLevel(); ++lev) {
      amrex::Print() << lev << ": " << amr2.Geom(lev).ProbDomain()
                     << amr2.Geom(lev).Domain() << amr2.Geom(lev).CellSize(0)
                     << "\n";
    }
    amr2.initialInit(0.0, 10000.0);
    

    // assert(dx1 = dx2)
    // assert(Domain lo = Domain lo)

    amrex::Print() << ">>> Filling in amr2 ...\n";
    int in_level = 0;
    Vector<int> offset_cell = {AMREX_D_DECL(0, 0, 0)};
    //ParmParse pp("manip");
    pp.query("in_level", in_level);
    // pp.queryarr("offset_cell", offset_cell, 0, AMREX_SPACEDIM); // TODO

    auto& amrlev1 = amr1.getLevel(in_level);
    auto& amrlev2 = amr2.getLevel(fill_level);
    MultiFab& mf1 = amrlev1.get_new_data(0); // State_Type
    MultiFab& mf2 = amrlev2.get_new_data(0);
    // amrlev2.state[0].allocOldData();
    // MultiFab& mfo1 = amrlev1.get_old_data(0);
    // MultiFab& mfo2 = amrlev2.get_old_data(0);
    MFIter::allowMultipleMFIters(true);
    for (MFIter mfi1(mf1); mfi1.isValid(); ++mfi1) {
      Box bx1 = mfi1.tilebox();
      for (MFIter mfi2(mf2); mfi2.isValid(); ++mfi2) {
        Box bx2 = mfi2.tilebox();
        Box bx = bx1 & bx2;

        if (bx.ok()) {
          auto const& arr1 = mf1.array(mfi1);
          auto const& arr2 = mf2.array(mfi2);
          // auto const& arro1 = mfo1.array(mfi1);
          // auto const& arro2 = mfo2.array(mfi2);
          ParallelFor(bx, LEN_STATE, [=](int i, int j, int k, int n) {
            arr2(i, j, k, n) = arr1(i, j, k, n);
            // arro2(i, j, k, n) = arro1(i, j, k, n);
          });
        }
      }
    }

    if (fill_level > 0) {
      for (int lev = fill_level - 1; lev >= 0; --lev)
        amr2.getLevel(lev).post_timestep(0);
    }

    amrex::Print() << ">>> Writing chkfile ...\n";
    amr2.checkPoint();

    amrex::Print() << ">>> Writing pltfile ...\n";
    amr2.writePlotFile();

    // #if CNS_USE_EB
    //     Geometry bigger_geom = (amr1.Geom(0).Domain() > amr2.Geom(0).Domain())
    //                              ? amr1.Geom(amr1.maxLevel())
    //                              : amr2.Geom(amr1.maxLevel());
    //     AmrLevel::SetEBSupportLevel(EBSupport::full);
    //     AmrLevel::SetEBMaxGrowCells(6, 6, 6);
    //     initialize_EB2(bigger_geom, amr1.maxLevel(), amr1.maxLevel());
    // #endif

    // amrex::Print() << ">>> Writing chkfile ...\n";
    // amr1.checkPoint();

    // amrex::Print() << ">>> Writing pltfile ...\n";
    // amr1.writePlotFile();

    // load_amr1();

    // init_amr2(start_time, stop_time);

    // fill_amr2();

    // output_amr2();
  }

  std::tuple<RealBox, int, Vector<int>> get_amr1()
  {
    RealBox rb_in;
    int max_level_in;
    Vector<int> n_cell_in(AMREX_SPACEDIM);

    ParmParse pp("amr");
    std::string filename;
    pp.get("restart", filename);
    std::string File(filename + "/Header");

    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);


//	std::string line;
//	for (int i = 0; i < 20; ++i) {
//		std::getline(is, line);
//		std::cout << i << ": " << line << '\n';
//	}

    int spdim;
    bool new_checkpoint_format = false;
    std::string first_line;
    std::getline(is, first_line);
    const std::string CheckPointVersion("CheckPointVersion_1.0");
    if (first_line == CheckPointVersion) {
      new_checkpoint_format = true;
      is >> spdim;
    } else {
      spdim = atoi(first_line.c_str());
    }
    // if (spdim != AMREX_SPACEDIM) {
    //   amrex::ErrorStream() << "Amr::restart(): bad spacedim = " << spdim << '\n';
    //   amrex::Abort();
    // }
    // amrex::Print() << "SpaceDim = " << spdim;

    Real cumtime;
    int max_lev;
    int finest_level;
    is >> cumtime;
    is >> max_lev;
    is >> finest_level;
    //amrex::Print() << ", cumtime = " << cumtime;
    //amrex::Print() << ", max_lev = " << max_lev;
    //amrex::Print() << ", finest_level = " << finest_level;

    Vector<Geometry> geom(max_lev + 1);
    //amrex::Print() << "\nGeom:\n";
    for (int lev = 0; lev <= max_lev; ++lev) {
      is  >> geom[lev];
      //amrex::Print() << geom[lev].ProbDomain() << geom[lev].Domain() << "\n";
    }

    rb_in = geom[0].ProbDomain();
    max_level_in = max_lev;
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      n_cell_in[i] = geom[0].Domain().bigEnd(i) - geom[0].Domain().smallEnd(i) + 1;
    }

    return {rb_in, max_level_in, n_cell_in};
  }

  std::tuple<RealBox, int, Vector<int>> get_amr2()
  {
    const int max_level_in = 0;

    ParmParse pp("amr");
    ParmParse ppg("geometry");
    Vector<Real> prob_lo, prob_hi;
    Vector<int> n_cell_in;
    ppg.getarr("prob_lo", prob_lo, 0, AMREX_SPACEDIM);
    ppg.getarr("prob_hi", prob_hi, 0, AMREX_SPACEDIM);
    pp.getarr("n_cell", n_cell_in, 0, AMREX_SPACEDIM);

    RealBox rb_in(prob_lo.data(), prob_hi.data());

    return {rb_in, max_level_in, n_cell_in};
  }

  // void load_amr1() {
  //   amr1.restart(amr1.restart_chkfile);
  // }

  // void init_amr2(Real strt_time, Real stop_time) {
  //   amr2.initialInit(strt_time, stop_time);
  //   amr2.updateInSitu();
  // }

  // void fill_amr2() {

  // }

  // void output_amr2() {
  //   amr2.checkPoint();
  //   amr2.writePlotFile();
  // }
};

int main(int argc, char* argv[])
{
  amrex::Initialize(argc, argv);

  TwoAmr driver;

  amrex::Finalize();
  return 0;
}
