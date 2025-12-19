#ifndef RIGIDBODYPROPERTIES_H_
#define RIGIDBODYPROPERTIES_H_

#include <AMReX.H>
#include <AMReX_RealVect.H>
#include <AMReX_Print.H>

// Include eib_cgal.h to get GeomType definition
#include <eib_cgal.h>

using namespace amrex;

namespace FSI {

struct RigidBodyProps {
    Real mass;
    RealVect xcenter;           // Center of Mass
    Real inertia[3][3];    // Inertia Tensor (3x3)
    Real volume;           // Volume (3D) or Area (2D)

    RigidBodyProps() {
        mass = 0.0;
        volume = 0.0;
        for(int i=0; i<AMREX_SPACEDIM; ++i) xcenter[i] = 0.0;
        for(int i=0; i<3; ++i)
            for(int j=0; j<3; ++j)
                inertia[i][j] = 0.0;
    }
};

class RigidBodyProperties {
public:

    // Compute properties for a given geometry and density
    static RigidBodyProps computeProperties(const GeomType& geom, Real rho) {
        RigidBodyProps props;

#if (AMREX_SPACEDIM == 2)
        // -------------------------------------------------------------------
        // 2D Implementation (Polygon)
        // -------------------------------------------------------------------
        // GeomType is Polygon2D (wrapper around CGAL::Polygon_2)
        // We can access vertices via geom.vertex(i)
        
        int n = geom.size();
        if (n < 3) return props;

        Real area = 0.0;
        Real cx = 0.0;
        Real cy = 0.0;
        Real Izz = 0.0;

        // Iterate over edges
        for (int i = 0; i < n; ++i) {
            const auto& p1 = geom.vertex(i);
            const auto& p2 = geom.vertex((i + 1) % n);

            Real x1 = CGAL::to_double(p1.x());
            Real y1 = CGAL::to_double(p1.y());
            Real x2 = CGAL::to_double(p2.x());
            Real y2 = CGAL::to_double(p2.y());

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

        // Izz about origin
        Izz /= 12.0;
        Izz = std::abs(Izz); // Ensure positive

        // Parallel axis theorem to move to CM
        // I_cm = I_origin - mass * d^2
        Real d2 = cx*cx + cy*cy;
        Real Izz_cm = (Izz * rho) - props.mass * d2; // Note: Izz above was geometric, multiply by rho

        // Fill inertia tensor (3x3 for compatibility)
        // In 2D, only Izz is relevant for rotation in plane
        props.inertia[2][2] = Izz_cm;

#elif (AMREX_SPACEDIM == 3)
        // -------------------------------------------------------------------
        // 3D Implementation (Polyhedron)
        // -------------------------------------------------------------------
        // GeomType is CGAL::Polyhedron_3
        // We iterate over facets and form tetrahedrons with the origin

        Real total_vol = 0.0;
        RealVect total_xcenter(0.0, 0.0, 0.0);
        
        // Covariance matrix of volume: C_ij = integral(x_i * x_j dV)
        Real C[3][3] = {{0.0}};

        // Iterate over facets
        for (auto f = geom.facets_begin(); f != geom.facets_end(); ++f) {
            // Get vertices of the triangle
            auto h = f->halfedge();
            const auto& p1_cgal = h->vertex()->point();
            const auto& p2_cgal = h->next()->vertex()->point();
            const auto& p3_cgal = h->next()->next()->vertex()->point();

            RealVect p1(CGAL::to_double(p1_cgal.x()), CGAL::to_double(p1_cgal.y()), CGAL::to_double(p1_cgal.z()));
            RealVect p2(CGAL::to_double(p2_cgal.x()), CGAL::to_double(p2_cgal.y()), CGAL::to_double(p2_cgal.z()));
            RealVect p3(CGAL::to_double(p3_cgal.x()), CGAL::to_double(p3_cgal.y()), CGAL::to_double(p3_cgal.z()));

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
};

} // namespace FSI

#endif // RIGIDBODYPROPERTIES_H_
