#ifndef RIGIDBODYPROPERTIES_H_
#define RIGIDBODYPROPERTIES_H_

#include <cmath>
#include <vector>
#include <AMReX.H>
#include <AMReX_RealVect.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

// GeomType is Polygon2D (2D) or TriMesh (3D), defined in ibm_bvh_defs.h
// Already included transitively via ibm_solver.h -> ibm_containers.h -> ibm_backend_bvh.h -> ibm_bvh.h -> ibm_bvh_defs.h
#include <ibm_bvh_defs.h>

using namespace amrex;

namespace FSI {

struct RigidBodyProps {
    Real mass;
    RealVect xcenter;           // Center of Mass
    Real inertia[3][3];         // Inertia Tensor (3x3, about CM)
    Real volume;                // Volume (3D) or Area (2D)
    Real rho_solid;             // Uniform solid density

    RigidBodyProps() {
        mass = 0.0;
        volume = 0.0;
        rho_solid = 0.0;
        for(int i=0; i<AMREX_SPACEDIM; ++i) xcenter[i] = 0.0;
        for(int i=0; i<3; ++i)
            for(int j=0; j<3; ++j)
                inertia[i][j] = 0.0;
    }
};

class RigidBodyProperties {
public:

    // Return rigid-body inertial properties for geometry geom_id.
    //
    // The result is computed once (on the first call) and cached; subsequent
    // calls return the cached value without re-reading inputs or re-integrating
    // the geometry mesh.  This is correct because mass, inertia, and the
    // reference centre of mass are constant properties of a rigid body.
    //
    // Two initialisation modes (selected via the fsi.* input namespace):
    //   Manual : fsi.mass > 0 is provided; remaining properties are read from
    //            the input file; the geometry mesh is not consulted.
    //   Auto   : fsi.mass absent; fsi.rho_solid (uniform density) required;
    //            mass, centre of mass, and inertia tensor are computed by
    //            numerical volume integration over the discretised geometry.
    //
    // Input keys:
    //   fsi.rho_solid       = <Real>            (auto: uniform solid density)
    //   fsi.mass            = <Real>            (manual: total mass)
    //   fsi.xcenter         = <Real ...>        (manual: centre of mass)
    //   fsi.Izz             = <Real>            (2D manual: polar moment about CM)
    //   fsi.inertia_diag    = <Real Real Real>  (3D manual: Ixx Iyy Izz about CM)
    //   fsi.inertia_offdiag = <Real Real Real>  (3D manual: Ixy Ixz Iyz about CM)
    static RigidBodyProps readOrCompute(const GeomType& geom, int geom_id = 0)
    {
        // Cache: indexed by geom_id; valid flag guards first-call initialisation.
        static std::vector<RigidBodyProps> cache;
        static std::vector<bool>          valid;

        // Grow cache on first encounter of a new geom_id
        if (geom_id >= static_cast<int>(cache.size())) {
            cache.resize(geom_id + 1);
            valid.resize(geom_id + 1, false);
        }
        if (valid[geom_id]) return cache[geom_id];

        amrex::ParmParse pp("fsi");

        Real manual_mass = -1.0;
        pp.query("mass", manual_mass);

        if (manual_mass > 0.0) {
            cache[geom_id] = readManual(pp, manual_mass);
        } else {
            Real rho = 0.0;
            if (!pp.query("rho_solid", rho) || rho <= 0.0) {
                amrex::Abort("FSI: fsi.rho_solid is required when fsi.mass is not provided");
            }
            cache[geom_id] = computeFromGeometry(geom, rho);
            amrex::Print() << "[FSI] Geometry " << geom_id
                           << ": computed from mesh (rho_solid=" << rho << ")\n"
                           << "  volume=" << cache[geom_id].volume
                           << "  mass="   << cache[geom_id].mass << "\n"
                           << "  xcenter=(" << AMREX_D_TERM(cache[geom_id].xcenter[0],
                                               << "," << cache[geom_id].xcenter[1],
                                               << "," << cache[geom_id].xcenter[2]) << ")\n"
                           << "  Izz=" << cache[geom_id].inertia[2][2] << "\n";
        }

        valid[geom_id] = true;
        return cache[geom_id];
    }

private:

    // Populate RigidBodyProps from input-file values (manual mode).
    static RigidBodyProps readManual(amrex::ParmParse& pp, Real mass) {
        RigidBodyProps props;
        props.mass = mass;

        // Center of mass (optional, default origin)
        std::vector<Real> xc;
        if (pp.queryarr("xcenter", xc)) {
            for (int d = 0; d < AMREX_SPACEDIM && d < static_cast<int>(xc.size()); ++d) {
                props.xcenter[d] = xc[d];
            }
        }

        // Volume (optional)
        pp.query("volume", props.volume);

        // Density (back-compute if volume known)
        if (props.volume > 0.0) {
            props.rho_solid = props.mass / props.volume;
        }

#if (AMREX_SPACEDIM == 2)
        // 2D: only Izz matters
        pp.query("Izz", props.inertia[2][2]);
#else
        // 3D: diagonal + off-diagonal
        std::vector<Real> diag, offdiag;
        if (pp.queryarr("inertia_diag", diag) && diag.size() >= 3) {
            props.inertia[0][0] = diag[0];
            props.inertia[1][1] = diag[1];
            props.inertia[2][2] = diag[2];
        }
        if (pp.queryarr("inertia_offdiag", offdiag) && offdiag.size() >= 3) {
            props.inertia[0][1] = props.inertia[1][0] = offdiag[0];
            props.inertia[0][2] = props.inertia[2][0] = offdiag[1];
            props.inertia[1][2] = props.inertia[2][1] = offdiag[2];
        }
#endif

        amrex::Print() << "[FSI] Manual input: mass=" << props.mass
                       << "  xcenter=(" << AMREX_D_TERM(props.xcenter[0],
                                          << "," << props.xcenter[1],
                                          << "," << props.xcenter[2]) << ")\n";
        return props;
    }

public:

    // Compute mass, centre of mass, and inertia tensor from the discretised
    // geometry assuming uniform solid density rho.  Exposed publicly so that
    // callers can bypass the input-file path when needed.
    static RigidBodyProps computeFromGeometry(const GeomType& geom, Real rho) {
        RigidBodyProps props;
        props.rho_solid = rho;

#if (AMREX_SPACEDIM == 2)
        // -------------------------------------------------------------------
        // 2D: signed-area shoelace formula for a closed polygon.
        // Centroid and second moment of area about the origin are accumulated
        // edge-by-edge; the parallel axis theorem shifts the result to the CM.
        // GeomType = Polygon2D: vertex access via geom.vertex(i)[0/1].
        // -------------------------------------------------------------------
        int n = static_cast<int>(geom.size());
        if (n < 3) return props;

        Real area = 0.0;
        Real cx = 0.0;
        Real cy = 0.0;
        Real Izz = 0.0;  // second moment of area about the origin

        for (int i = 0; i < n; ++i) {
            const auto& p1 = geom.vertex(i);
            const auto& p2 = geom.vertex((i + 1) % n);

            Real x1 = p1[0];
            Real y1 = p1[1];
            Real x2 = p2[0];
            Real y2 = p2[1];

            Real cross = x1 * y2 - x2 * y1;
            area += cross;

            cx += (x1 + x2) * cross;
            cy += (y1 + y2) * cross;

            // Polar second moment about the origin:
            // integral(x^2+y^2)dA = sum_edges cross*(x1^2+x1*x2+x2^2+y1^2+y1*y2+y2^2)/12
            Izz += cross * (x1*x1 + x1*x2 + x2*x2 + y1*y1 + y1*y2 + y2*y2);
        }

        area *= 0.5;
        props.volume = std::abs(area);  // area in 2D
        props.mass   = props.volume * rho;

        if (std::abs(area) > 1.0e-12) {
            cx /= (6.0 * area);
            cy /= (6.0 * area);
        }
        props.xcenter[0] = cx;
        props.xcenter[1] = cy;

        // Shift second moment to CM via the parallel axis theorem:
        //   I_cm = rho * I_geo_origin - mass * |x_cm|^2
        Izz = std::abs(Izz) / 12.0;
        props.inertia[2][2] = rho * Izz - props.mass * (cx*cx + cy*cy);

#elif (AMREX_SPACEDIM == 3)
        // -------------------------------------------------------------------
        // 3D: signed-volume decomposition into tetrahedra with the origin as apex.
        // Volume, centroid, and the symmetric covariance matrix C_ij are
        // accumulated face-by-face; the full inertia tensor about the origin is
        // assembled from C, then shifted to the CM via the parallel axis theorem.
        // GeomType = TriMesh: faces[fi][0..2] are vertex indices into vertices[].
        // -------------------------------------------------------------------

        int nf = geom.num_faces();
        if (nf < 1) return props;

        Real total_vol = 0.0;
        RealVect total_xcenter(0.0, 0.0, 0.0);
        
        // Covariance matrix of volume: C_ij = integral(x_i * x_j dV)
        Real C[3][3] = {{0.0}};

        // Iterate over triangular faces
        for (int fi = 0; fi < nf; ++fi) {
            const auto& face = geom.faces[fi];
            const auto& v0 = geom.vertices[face[0]];
            const auto& v1 = geom.vertices[face[1]];
            const auto& v2 = geom.vertices[face[2]];

            RealVect p1(v0[0], v0[1], v0[2]);
            RealVect p2(v1[0], v1[1], v1[2]);
            RealVect p3(v2[0], v2[1], v2[2]);

            // Signed volume of tetrahedron O-P1-P2-P3
            // V = 1/6 * det(p1, p2, p3)
            Real det = p1[0]*(p2[1]*p3[2] - p2[2]*p3[1]) 
                     - p1[1]*(p2[0]*p3[2] - p2[2]*p3[0]) 
                     + p1[2]*(p2[0]*p3[1] - p2[1]*p3[0]);
            Real vol = det / 6.0;

            total_vol += vol;

            // Centroid of tetrahedron
            RealVect tet_xcenter = (p1 + p2 + p3) * 0.25; // Origin is (0,0,0)
            total_xcenter += tet_xcenter * vol;

            // Covariance integrals for the tetrahedron (Dobrovolskis 1996):
            //   C_ii += V/10  * (pi^2 + pi*pj + pj^2)  [diagonal, i==j]
            //   C_ij += V/20  * (2pi_a*pj_a + cross terms) [off-diagonal]
            
            for(int i=0; i<3; ++i) {
                for(int j=i; j<3; ++j) {
                    Real term = 0.0;
                    if (i == j) {
                        term = (p1[i]*p1[i] + p2[i]*p2[i] + p3[i]*p3[i] + 
                                p1[i]*p2[i] + p2[i]*p3[i] + p3[i]*p1[i]);
                        term *= (vol / 10.0);
                    } else {
                        term = (2*p1[i]*p1[j] + p1[i]*p2[j] + p1[i]*p3[j] +
                                p2[i]*p1[j] + 2*p2[i]*p2[j] + p2[i]*p3[j] +
                                p3[i]*p1[j] + p3[i]*p2[j] + 2*p3[i]*p3[j]);
                        term *= (vol / 20.0);
                    }
                    C[i][j] += term;
                }
            }
        }

        // Fill symmetric part of C
        C[1][0] = C[0][1]; C[2][0] = C[0][2]; C[2][1] = C[1][2];

        props.volume = std::abs(total_vol);
        props.mass = props.volume * rho;

        if (std::abs(total_vol) > 1e-12) {
            total_xcenter /= total_vol;
        }
        props.xcenter = total_xcenter;

        // Inertia tensor about the origin: I_ij = rho*(delta_ij*tr(C) - C_ij)
        Real I_origin[3][3];
        I_origin[0][0] = rho * (C[1][1] + C[2][2]);
        I_origin[1][1] = rho * (C[0][0] + C[2][2]);
        I_origin[2][2] = rho * (C[0][0] + C[1][1]);
        I_origin[0][1] = I_origin[1][0] = -rho * C[0][1];
        I_origin[0][2] = I_origin[2][0] = -rho * C[0][2];
        I_origin[1][2] = I_origin[2][1] = -rho * C[1][2];

        // Parallel axis theorem: I_cm = I_origin - M*(|d|^2 δ - d⊗d)
        const Real cx = props.xcenter[0];
        const Real cy = props.xcenter[1];
        const Real cz = props.xcenter[2];

        props.inertia[0][0] = I_origin[0][0] - props.mass * (cy*cy + cz*cz);
        props.inertia[1][1] = I_origin[1][1] - props.mass * (cx*cx + cz*cz);
        props.inertia[2][2] = I_origin[2][2] - props.mass * (cx*cx + cy*cy);
        props.inertia[0][1] = props.inertia[1][0] = I_origin[0][1] + props.mass * cx * cy;
        props.inertia[0][2] = props.inertia[2][0] = I_origin[0][2] + props.mass * cx * cz;
        props.inertia[1][2] = props.inertia[2][1] = I_origin[1][2] + props.mass * cy * cz;

#endif

        return props;
    }

}; // class RigidBodyProperties

} // namespace FSI

#endif // RIGIDBODYPROPERTIES_H_
