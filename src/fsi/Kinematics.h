#ifndef KINEMATICS_H_
#define KINEMATICS_H_

#include <AMReX.H>
#include <AMReX_RealVect.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

// Include CNS to access IBM::ib
#include <CNS.h> 

using namespace amrex;

namespace FSI {

class Kinematics {
public:
    struct ForceMoment {
        RealVect force;
        RealVect moment; 
    };

    // Calculate loads for a specific geometry
    // geom_id: Index of the geometry to calculate loads for
    // ref_point: Reference point for moment calculation
    static ForceMoment computeLoads(int geom_id, const RealVect& ref_point) {
        
        ForceMoment loads;
        loads.force = RealVect(AMREX_D_DECL(0.0, 0.0, 0.0));
        loads.moment = RealVect(AMREX_D_DECL(0.0, 0.0, 0.0));

        // Access the global IBM instance
        auto& ib = IBM::ib;
        
        // Check if geometry index is valid
        if (geom_id < 0 || geom_id >= ib.ngeom) {
            amrex::Print() << "Error: Invalid geometry index " << geom_id << " in computeLoads.\n";
            return loads;
        }

        // Get the range of faces for this geometry
        int start_idx = ib.geom_offsets[geom_id];
        int end_idx   = ib.geom_offsets[geom_id + 1];

        // Loop over faces
        for (int i = start_idx; i < end_idx; ++i) {
            
            // 1. Check ownership: Only calculate for faces owned by this rank
            // elemfound: 1 = owned/valid
            if (ib.surfphys_soa.elemfound[i] != 1) continue;

            // 2. Retrieve Physical Data
            Real p    = ib.surfphys_soa.pressure[i];
            Real tau1 = ib.surfphys_soa.tau1[i];
            Real tau2 = 0.0;
#if (AMREX_SPACEDIM == 3)
            tau2 = ib.surfphys_soa.tau2[i];
#endif

            // 3. Retrieve Geometric Data
            Real area = ib.SurfElem_a[i].measure;
            
            RealVect centroid;
            for(int d=0; d<AMREX_SPACEDIM; ++d) {
                centroid[d] = ib.SurfElem_a[i].centroid[d];
            }

            RealVect normal;
            RealVect tangent1;
#if (AMREX_SPACEDIM == 3)
            RealVect tangent2;
#endif
            for(int d=0; d<AMREX_SPACEDIM; ++d) {
                normal[d]   = ib.LocalFrame_a[i].normal[d];
                tangent1[d] = ib.LocalFrame_a[i].tangent1[d];
#if (AMREX_SPACEDIM == 3)
                tangent2[d] = ib.LocalFrame_a[i].tangent2[d];
#endif
            }

            // 4. Calculate Elemental Force dF
            // Pressure acts into the surface (-n)
            // Shear stress acts along tangents
            RealVect dF = -1.0 * p * normal * area; 
            dF += tau1 * tangent1 * area;
#if (AMREX_SPACEDIM == 3)
            dF += tau2 * tangent2 * area;
#endif

            // 5. Calculate Elemental Moment dT = r x dF
            RealVect r = centroid - ref_point;
            
            // Cross Product Logic
#if (AMREX_SPACEDIM == 3)
            RealVect dT = r.crossProduct(dF);
            loads.moment += dT;
#else
            // 2D Cross Product (Scalar result in Z direction)
            // r = (rx, ry), dF = (fx, fy)
            // Mz = rx * fy - ry * fx
            Real Mz = r[0] * dF[1] - r[1] * dF[0];
            // Store in first component for 2D
            loads.moment[0] += Mz; 
#endif
            loads.force += dF;
        }

        // 6. MPI Reduction
        ParallelDescriptor::ReduceRealSum(loads.force.begin(), AMREX_SPACEDIM);
        
#if (AMREX_SPACEDIM == 3)
        ParallelDescriptor::ReduceRealSum(loads.moment.begin(), AMREX_SPACEDIM);
#else
        // In 2D, we only used moment[0]
        ParallelDescriptor::ReduceRealSum(loads.moment.begin(), 1);
#endif

        return loads;
    }
};

} // namespace FSI

#endif // KINEMATICS_H_
