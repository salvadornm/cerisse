#ifndef RIGIDBODYPROPERTIES_H_
#define RIGIDBODYPROPERTIES_H_

#include <AMReX.H>
#include <AMReX_RealVect.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

// GeomType is Polygon2D (2D) or TriMesh (3D), defined in bvh_types.h
// Already included transitively via eib.h -> eib_data.h -> eib_bvh.h -> bvh.h -> bvh_types.h
#include <bvh_types.h>

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

    // ----------------------------------------------------------------
    // Read properties from input file (fsi.* namespace).
    // If fsi.mass is provided, all values come from inputs (skip geometry
    // computation).  Otherwise compute from geometry assuming uniform
    // density fsi.rho_solid (required).
    //
    // Input keys (per-geometry, indexed by geom_id):
    //   fsi.rho_solid       = <Real>          (uniform density, required for auto)
    //   fsi.mass            = <Real>          (manual override — triggers skip)
    //   fsi.xcenter         = <Real Real ...> (manual CM)
    //   fsi.Izz             = <Real>          (2D: polar moment about CM)
    //   fsi.inertia_diag    = <Real Real Real> (3D: Ixx Iyy Izz about CM)
    //   fsi.inertia_offdiag = <Real Real Real> (3D: Ixy Ixz Iyz about CM)
    // ----------------------------------------------------------------
    static RigidBodyProps readOrCompute(const GeomType& geom, int geom_id = 0) {

        amrex::ParmParse pp("fsi");

        // --- Check for manual override ---
        Real manual_mass = -1.0;
        pp.query("mass", manual_mass);

        if (manual_mass > 0.0) {
            return readManual(pp, manual_mass);
        }

        // --- Auto-compute from geometry ---
        Real rho = 0.0;
        if (!pp.query("rho_solid", rho) || rho <= 0.0) {
            amrex::Abort("FSI: fsi.rho_solid is required when fsi.mass is not provided");
        }

        RigidBodyProps props = computeFromGeometry(geom, rho);

        amrex::Print() << "[FSI] Geometry " << geom_id
                       << ": computed from mesh (rho_solid=" << rho << ")\n"
                       << "  volume=" << props.volume
                       << "  mass=" << props.mass << "\n"
                       << "  xcenter=(" << AMREX_D_TERM(props.xcenter[0],
                                          << "," << props.xcenter[1],
                                          << "," << props.xcenter[2]) << ")\n"
                       << "  Izz=" << props.inertia[2][2] << "\n";
        return props;
    }

private:

    // Read all properties from input file (manual mode)
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

    // Compute properties from geometry assuming uniform density rho.
    // Public so callers can still use it directly if needed.
    static RigidBodyProps computeFromGeometry(const GeomType& geom, Real rho) {
        RigidBodyProps props;
        props.rho_solid = rho;

#if (AMREX_SPACEDIM == 2)
        // -------------------------------------------------------------------
        // 2D Implementation (Polygon)
        // -------------------------------------------------------------------
        // GeomType is Polygon2D (bvh_types.h): vertices as GpuArray<Real,2>
        // Access via geom.vertex(i)[0], geom.vertex(i)[1]
        
        int n = static_cast<int>(geom.size());
        if (n < 3) return props;

        Real area = 0.0;
        Real cx = 0.0;
        Real cy = 0.0;
        Real Izz = 0.0;

        // Iterate over edges (shoelace formula)
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

            // Contribution to polar moment of inertia about origin
            // I_origin = integral (x^2 + y^2) dA
            // Formula: (x1*y2 - x2*y1) * (x1^2 + x1*x2 + x2^2 + y1^2 + y1*y2 + y2^2) / 12
            Izz += cross * (x1*x1 + x1*x2 + x2*x2 + y1*y1 + y1*y2 + y2*y2);
        }

        area *= 0.5;
        props.volume = std::abs(area); // Area in 2D
        props.mass = props.volume * rho;

        if (std::abs(area) > 1e-12) {
            cx /= (6.0 * area);
            cy /= (6.0 * area);
        }
        props.xcenter[0] = cx;
        props.xcenter[1] = cy;

        // Izz about origin (geometric second moment, needs rho to get inertia)
        Izz /= 12.0;
        Izz = std::abs(Izz);

        // Parallel axis theorem: I_cm = I_origin - mass * d^2
        Real d2 = cx*cx + cy*cy;
        Real Izz_cm = (Izz * rho) - props.mass * d2;

        // In 2D, only Izz is relevant for rotation in plane
        props.inertia[2][2] = Izz_cm;

#elif (AMREX_SPACEDIM == 3)
        // -------------------------------------------------------------------
        // 3D Implementation (TriMesh)
        // -------------------------------------------------------------------
        // GeomType is TriMesh (bvh_types.h): vertices[] and faces[] arrays
        // Each face has 3 vertex indices. Form tetrahedra with origin.

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

            // Inertia integrals (Covariance terms)
            // C_xx = integral(x^2 dV)
            // Formula for tetrahedron at origin:
            // integral(x^2) = V/10 * (x1^2 + x2^2 + x3^2 + x1x2 + x2x3 + x3x1)
            // integral(xy)  = V/20 * (2x1y1 + x1y2 + x1y3 + x2y1 + 2x2y2 + x2y3 + x3y1 + x3y2 + 2x3y3)
            
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

        // Compute Inertia Tensor about Origin
        // I_xx = integral(y^2 + z^2) = C_yy + C_zz
        // I_xy = -integral(xy) = -C_xy
        Real I_origin[3][3];
        I_origin[0][0] = rho * (C[1][1] + C[2][2]);
        I_origin[1][1] = rho * (C[0][0] + C[2][2]);
        I_origin[2][2] = rho * (C[0][0] + C[1][1]);
        
        I_origin[0][1] = I_origin[1][0] = -rho * C[0][1];
        I_origin[0][2] = I_origin[2][0] = -rho * C[0][2];
        I_origin[1][2] = I_origin[2][1] = -rho * C[1][2];

        // Parallel Axis Theorem to move to CM
        // I_cm = I_origin - M * (d^2 I - d tensor d)
        // I_cm_xx = I_origin_xx - M * (dy^2 + dz^2)
        // I_cm_xy = I_origin_xy - M * (-dx * dy) = I_origin_xy + M * dx * dy
        
        Real cx = props.xcenter[0];
        Real cy = props.xcenter[1];
        Real cz = props.xcenter[2];

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
