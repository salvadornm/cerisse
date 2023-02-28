#include <IBM.H>
#include <CGAL/Side_of_triangle_mesh.h>

using namespace amrex;
using namespace IBM;

IBFab::IBFab (const Box& b, int ncomp, bool alloc, bool shared, Arena* ar)    
              : BaseFab<bool>(b,ncomp,alloc,shared,ar) {}
IBFab::IBFab (const IBFab& rhs, MakeType make_type, int scomp, int ncomp) 
              : BaseFab<bool>(rhs,make_type,scomp,ncomp) {}
IBFab::~IBFab () { delete[] gpArray; }

void IBFab::allocateGPs(int numGPs) {
  delete[] gpArray;
  gpArray = new gp[numGPs]; 
  ngps = numGPs; }

IBMultiFab::IBMultiFab ( const BoxArray& bxs, const DistributionMapping& dm, 
                        const int nvar, const int ngrow, const MFInfo& info, 
                        const FabFactory<IBFab>& factory )  :
                        FabArray<IBFab>(bxs,dm,nvar,ngrow,info,factory) {}
IBMultiFab::~IBMultiFab () {}

// constructor without geometry init
// IB::IB (const Vector<BoxArray>& bxs, 
//         const Vector<DistributionMapping>& dm, 
//         const int nvar, const int nghost, const int max_level) {
IB::IB (Amr* pointer_amr, const int nvar, const int nghost, const int max_level, const Vector<GpuArray<Real,AMREX_SPACEDIM>> dx) {
        
        // store pointer to main Amr class object's instance
        IB::pamr = pointer_amr ;

        // Since AmrInfo class is protected! I believe this prevents access to its members here, like max_level,etc. So store basic AmrInfo for later easy access.
        max_lev = max_level;
        ref_ratio = pamr->refRatio();
        cellSizes.resize(max_level+1);
        cellSizes = dx;

        // create IBMultiFabs at each level and store pointers to it
        mfa->resize(max_level+1);
        for (int lev=0; lev<=max_level; lev++) {
          mfa->at(lev) = new IBMultiFab(pamr->boxArray(lev),pamr->DistributionMap(lev),nvar,nghost);
        } }

IB::~IB () { 
  // For a vector of pointers, we must delete each pointer to object individually. https://shorturl.at/dHPXZ 
 for (int lev=0; lev<=mfa->size(); lev++) {
          delete[] mfa->at(lev);
        }
}

void IB::compute_markers () {

  // int nb_inside = 0;
  // int nb_boundary = 0;
  // CGAL::Side_of_triangle_mesh<Polyhedron, K> inside(IB::geom);
  // // std::vector<Point> points;
  // for (std::size_t i = 0; i < nb_points; ++i)
  // {
  //   CGAL::Bounded_side res = inside(points[i]);
  //   if (res == CGAL::ON_BOUNDED_SIDE) { ++nb_inside; }
  //   if (res == CGAL::ON_BOUNDARY) { ++nb_boundary; }
  // }

  // std::cerr << "Total query size: " << points.size() << std::endl;
  // std::cerr << "  " << nb_inside << " points inside " << std::endl;
  // std::cerr << "  " << nb_boundary << " points on boundary " << std::endl;
  // std::cerr << "  " << points.size() - nb_inside - nb_boundary << " points outside " << std::endl;

        // pamr.max_grid_size()

  for (int lev=0; lev<mfa->size(); lev++) {
    IBMultiFab *mfab = mfa->at(lev);
    for (MFIter mfi(*mfab,false); mfi.isValid(); ++mfi) {
        // IBM::IBFab &fab = (*mfab).[mfi];
        IBM::IBFab &fab = mfab->get(mfi);
        const int *lo = fab.loVect();
        const int *hi = fab.hiVect();
        amrex::Print() << "Current level/Max level = " << lev  << "/" << mfa->size()-1 << std::endl;
        amrex::Print() << lo[0] << " " << lo[1] << " " << lo[2] << std::endl;
        amrex::Print() << hi[0] << " " << hi[1] << " " << hi[2] << std::endl;


        fab.setVal(false); // initialise
        const Box& bx = mfi.fabbox(); // box with ghost points
        auto const bool_array = fab.array(); // boolean array

        // amrex::ParallelFor(bx,
        // [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        // {
        //   amrex::Print() << i << " "<< j << " " << k << " " << std::endl;
        //   amrex::Print() << bool_array(i,j,k,1) << " " << std::endl;
        //   fab.
        // });
        
        // const Box& bx = mfi.validbox(); // box without ghost points
        // auto const& sfab = (*mfab).array(mfi);
          // bool_array.begin
          // bool_array.ncomp 


        amrex::Print() << "fab ngps = " << fab.ngps << std::endl;
        amrex::Print() << "------------------- " << std::endl;
    }
  }
}


void IB::read_geom(const std::string filename) {
  namespace PMP = CGAL::Polygon_mesh_processing;

  if(!PMP::IO::read_polygon_mesh(filename, IB::geom) || CGAL::is_empty(IB::geom) || !CGAL::is_triangle_mesh(IB::geom))
  {
    std::cerr << "Invalid input." << std::endl;
    exit(1);
  }
    amrex::Print() << "----------------------------------" << std::endl;
    amrex::Print() << "Geometry " << filename << " read"<< std::endl;
    amrex::Print() << "----------------------------------" << std::endl;
}
























///////////////////////////////CHECK IBMULTIFABS///////////////////
// int temp=0 ;
// for (MFIter mfi(*mfa[0],false); mfi.isValid(); ++mfi) // without tiling
// {
//     IBM::IBFab &fab = (*mfa[0])[mfi];
//     const int *lo = fab.loVect();
//     const int* hi = fab.hiVect();

//     amrex::Print() << "Level " << max_level << std::endl;
//     amrex::Print() << lo[0] << " " << lo[1] << " " << lo[2] << std::endl;
//     amrex::Print() << hi[0] << " " << hi[1] << " " << hi[2] << std::endl;

    
//     amrex::Print() << "fab allocated " << fab.isAllocated() << std::endl;
//     fab.allocateGPs(temp);
//     amrex::Print() << "fab ngps = " << fab.ngps << std::endl;

//     amrex::Print() << "------------------- " << std::endl;
//     temp += 1;
// }

///////////////////////////////////////////NORMAL MultiFab//////////////////////
  // MultiFab data( get_new_data(State_Type).boxArray(), get_new_data(State_Type).DistributionMap(), NUM_STATE, 1, MFInfo(), Factory());

  // const amrex::Real cur_time = state[State_Type].curTime();
  // FillPatch(*this, data, data.nGrow(), cur_time, State_Type, Density, NUM_STATE, 0);

  //   for (MFIter mfi(data,TilingIfNotGPU()); mfi.isValid(); ++mfi)
  //   {
  //       FArrayBox& fab = data[mfi];
  //       const int *lo = fab.loVect();
  //       const int* hi = fab.hiVect();

  //       amrex::Print() << "Level " << level << std::endl;
  //       amrex::Print() << lo[0] << " " << lo[1] << " " << lo[2] << std::endl;
  //       amrex::Print() << hi[0] << " " << hi[1] << " " << hi[2] << std::endl;
        
  //   }
////////////////////////////////////////////////////////////////////////////////


// 2. compute and save solid and ghost point field


// 3. compute ghost point state


// 4. compute surface data