#ifndef EIB_IO_H_
#define EIB_IO_H_

// ============================================================================
// eib_io.h — Geometry I/O, VTK output, MPI gather, and CSR builder
//
// This file contains in-class definitions for member functions of eib_t.
// It is #include'd INSIDE the class body in eib.h and must NOT be included
// independently.
//
// Contents:
//   1. gatherSurfData   — MPI Gatherv of surface physical data to rank 0
//   2. plotSURF          — Write VTK PolyData (.vtp) files per geometry
//   3. buildCSR          — Compressed Sparse Row structure for face iteration
//   4. getGeomIdx        — Map global element index to geometry index
//   5. read_geom         — Load geometry files, build BVH, fill geometry cache
// ============================================================================

// ============================================================================
// 1. gatherSurfData
// ============================================================================

public:

void gatherSurfData()
{
    int nprocs = amrex::ParallelDescriptor::NProcs();
    if (nprocs == 1) return; // Serial execution: data is already on Rank 0

    int my_rank = amrex::ParallelDescriptor::MyProc();
    MPI_Comm comm = amrex::ParallelDescriptor::Communicator();

    if (ntotalfaces <= 0) {
        amrex::Abort("gatherSurfData: ntotalfaces must be positive");
    }

    // Structure for packing data to send to Rank 0
    // This struct must be trivially copyable to be safely sent via MPI as raw bytes.
    struct SurfOut {
        int elem;      // face index (global)
        int lev;       // AMR level
        int rank;      // Owning rank
        int ipq;       // Interpolation quality (number of fluid points)
        Real p, T, dTdn, tau1, tau2; // Physical quantities
    };
    static_assert(std::is_trivially_copyable<SurfOut>::value,
                  "SurfOut must be trivially copyable for MPI");

    // ======================================================================
    // Step 1: Pack local owned faces
    // ======================================================================
    std::vector<SurfOut> send;
    send.reserve(surfphys_soa.filled_elems > 0 ? surfphys_soa.filled_elems : 1024);

    for (int i = 0; i < ntotalfaces; ++i) {
        if (surfphys_soa.elemfound[i] != 1) continue;

        SurfOut s{};
        s.elem = surfphys_soa.elemIdx[i];
        s.lev  = surfphys_soa.lev[i];
        s.rank = surfphys_soa.rank[i];
        s.ipq  = surfphys_soa.ip_quality[i];
        s.p    = surfphys_soa.pressure[i];
        s.T    = surfphys_soa.temperature[i];
        s.dTdn = surfphys_soa.dTdn[i];
        s.tau1 = surfphys_soa.tau1[i];
        s.tau2 = surfphys_soa.tau2[i];
        send.push_back(s);
    }

    int nlocal = static_cast<int>(send.size());

    // ======================================================================
    // Step 2: Gather counts on root
    // ======================================================================
    std::vector<int> counts = amrex::ParallelDescriptor::Gather(nlocal, 0);

    // ======================================================================
    // Step 3: Compute displacements & total size
    // ======================================================================
    std::vector<int> displs;
    int total = 0;
    if (my_rank == 0) {
        displs.resize(nprocs, 0);
        for (int r = 1; r < nprocs; ++r) displs[r] = displs[r - 1] + counts[r - 1];
        for (int r = 0; r < nprocs; ++r) total += counts[r];
    }

    // ======================================================================
    // Step 4: Gatherv bytes (robust, no custom MPI datatype)
    // ======================================================================
    const int typesize = sizeof(SurfOut);

    std::vector<char> sendbuf(reinterpret_cast<char*>(send.data()),
                              reinterpret_cast<char*>(send.data()) + nlocal * typesize);

    std::vector<char> recvbuf;
    std::vector<int> counts_b, displs_b;
    char* recvptr = nullptr;

    if (my_rank == 0) {
        recvbuf.resize(total * typesize);
        recvptr = recvbuf.data();

        counts_b.resize(nprocs);
        displs_b.resize(nprocs);
        for (int r = 0; r < nprocs; ++r) {
            if (counts[r] < 0) amrex::Abort("gatherSurfData: negative count detected");
            if (counts[r] > std::numeric_limits<int>::max() / typesize) {
                amrex::Abort("gatherSurfData: Gatherv byte count overflow");
            }
            counts_b[r] = counts[r] * typesize;
            displs_b[r] = displs[r] * typesize;
        }
    }

    amrex::ParallelDescriptor::Gatherv(sendbuf.data(), nlocal * typesize,
                                       recvptr, counts_b, displs_b, 0);

    // ======================================================================
    // Step 5: Unpack on root back into SoA
    // ======================================================================
    if (my_rank == 0) {
        auto* rec = reinterpret_cast<const SurfOut*>(recvbuf.data());

        for (int k = 0; k < total; ++k) {
            const auto& s = rec[k];
            const int i = s.elem;

            if (i < 0 || i >= ntotalfaces) {
                amrex::Abort("gatherSurfData: received face index out of range");
            }

            surfphys_soa.lev[i]        = s.lev;
            surfphys_soa.rank[i]       = s.rank;
            surfphys_soa.ip_quality[i] = s.ipq;

            surfphys_soa.pressure[i]    = s.p;
            surfphys_soa.temperature[i] = s.T;
            surfphys_soa.dTdn[i]        = s.dTdn;
            surfphys_soa.tau1[i]        = s.tau1;
            surfphys_soa.tau2[i]        = s.tau2;
        }
    }
}

// ============================================================================
// 2. plotSURF
// ============================================================================

void plotSURF(
    const amrex::Real time, int step, const std::string& prefix)
{
    // Only Rank 0 writes the file (serial I/O)
    if (amrex::ParallelDescriptor::MyProc() != 0) return;

    // Parse prefix to get base directory and filename prefix
    std::string base_dir = ".";
    std::string file_prefix = prefix;

    auto pos = prefix.find_last_of("/\\");
    if (pos != std::string::npos) {
        base_dir = prefix.substr(0, pos);
        file_prefix = prefix.substr(pos + 1);

        if (!amrex::UtilCreateDirectory(base_dir, 0755)) {
            amrex::Print() << "Error: Could not create directory " << base_dir << "\n";
        }
    }

    for (int i = 0; i < ngeom; ++i) {

        // Use geometry name from input file if available, fallback to geom{i}
        std::string gname = (i < static_cast<int>(geom_names.size()) && !geom_names[i].empty())
                            ? geom_names[i] : "geom" + std::to_string(i);

        // Create geometry specific directory: base_dir/{gname}
        std::string geom_dir = base_dir + "/" + gname;
        if (!amrex::UtilCreateDirectory(geom_dir, 0755)) {
            amrex::Print() << "Error: Could not create directory " << geom_dir << "\n";
        }

        // Construct filename: geom_dir/{gname}_{step}.vtp
        std::string filename = amrex::Concatenate(
            geom_dir + "/" + gname + "_", step, 5) + ".vtp";
        amrex::Print() << "Writing surface data for geometry " << i << " to " << filename << " ...\n";

        std::ofstream ofs(filename);
        if (!ofs.good()) {
            amrex::Print() << "Error: Cannot open file " << filename << " for writing.\n";
            continue;
        }

        // ==================================================================
        // 1. Count Points and Cells
        // ==================================================================
        long long n_points = 0;
        long long n_cells  = 0;
#if (AMREX_SPACEDIM == 3)
        n_points = geom_a[i].size_of_vertices();
        n_cells  = geom_a[i].size_of_facets();
#else
        n_points = geom_a[i].size();
        n_cells  = geom_a[i].size();
#endif

        // ==================================================================
        // 2. Write VTK XML Header
        // ==================================================================
        ofs << "<?xml version=\"1.0\"?>\n";
        ofs << "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
        ofs << "  <PolyData>\n";
        ofs << "    <Piece NumberOfPoints=\"" << n_points << "\" NumberOfPolys=\"" << n_cells << "\">\n";

        // ==================================================================
        // 3. Write Points (Vertices)
        // ==================================================================
        ofs << "      <Points>\n";
        ofs << "        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n";
#if (AMREX_SPACEDIM == 3)
#ifdef AMREX_USE_CGAL
        for (auto vit = geom_a[i].vertices_begin(); vit != geom_a[i].vertices_end(); ++vit) {
            const auto& p = vit->point();
            ofs << p.x() << " " << p.y() << " " << p.z() << " ";
        }
#else
        for (int v = 0; v < geom_a[i].num_vertices(); ++v) {
            ofs << geom_a[i].vertices[v][0] << " "
                << geom_a[i].vertices[v][1] << " "
                << geom_a[i].vertices[v][2] << " ";
        }
#endif
#else
        for (auto v = geom_a[i].vertices_begin(); v != geom_a[i].vertices_end(); ++v) {
            ofs << (*v)[0] << " " << (*v)[1] << " 0.0 ";
        }
#endif
        ofs << "\n        </DataArray>\n";
        ofs << "      </Points>\n";

        // ==================================================================
        // 4. Write Polys (Connectivity and Offsets)
        // ==================================================================
        ofs << "      <Polys>\n";

        ofs << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
#if (AMREX_SPACEDIM == 3)
#ifdef AMREX_USE_CGAL
        {
            // Build vertex index map for CGAL Polyhedron_3
            std::map<const void*, int> vindex;
            int idx = 0;
            for (auto vit = geom_a[i].vertices_begin(); vit != geom_a[i].vertices_end(); ++vit) {
                vindex[&(*vit)] = idx++;
            }
            for (auto fit = geom_a[i].facets_begin(); fit != geom_a[i].facets_end(); ++fit) {
                auto h = fit->facet_begin();
                ofs << vindex[&(*h->vertex())] << " ";
                ++h;
                ofs << vindex[&(*h->vertex())] << " ";
                ++h;
                ofs << vindex[&(*h->vertex())] << " ";
            }
        }
#else
        for (int f = 0; f < geom_a[i].num_faces(); ++f) {
            ofs << geom_a[i].faces[f][0] << " "
                << geom_a[i].faces[f][1] << " "
                << geom_a[i].faces[f][2] << " ";
        }
#endif
#else
        for (int k = 0; k < n_points; ++k) {
            int v1 = k;
            int v2 = (k + 1) % n_points;
            ofs << v1 << " " << v2 << " ";
        }
#endif
        ofs << "\n        </DataArray>\n";

        ofs << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
        long long current_offset = 0;
#if (AMREX_SPACEDIM == 3)
        for (int f = 0; f < n_cells; ++f) {
            current_offset += 3;
            ofs << current_offset << " ";
        }
#else
        for (int k = 0; k < n_cells; ++k) {
            current_offset += 2;
            ofs << current_offset << " ";
        }
#endif
        ofs << "\n        </DataArray>\n";
        ofs << "      </Polys>\n";

        // ==================================================================
        // 5. Write Cell Data (Physical Properties)
        // ==================================================================
        ofs << "      <CellData Scalars=\"Pressure\">\n";

        auto write_scalar_field = [&](const std::string& name, const amrex::Gpu::ManagedVector<Real>& field) {
            ofs << "        <DataArray type=\"Float64\" Name=\"" << name << "\" format=\"ascii\">\n";
            int start_idx = geom_offsets[i];
            if (start_idx + n_cells > static_cast<long long>(field.size())) {
                amrex::Print() << "Warning: Data field " << name << " size mismatch. Writing zeros.\n";
                for (int k = 0; k < n_cells; ++k) ofs << "0.0 ";
            } else {
                for (int k = 0; k < n_cells; ++k) {
                    ofs << field[start_idx + k] << " ";
                }
            }
            ofs << "\n        </DataArray>\n";
        };

        auto write_int_field = [&](const std::string& name, const amrex::Gpu::ManagedVector<int>& field) {
            ofs << "        <DataArray type=\"Int32\" Name=\"" << name << "\" format=\"ascii\">\n";
            int start_idx = geom_offsets[i];
            if (start_idx + n_cells > static_cast<long long>(field.size())) {
                amrex::Print() << "Warning: Data field " << name << " size mismatch. Writing zeros.\n";
                for (int k = 0; k < n_cells; ++k) ofs << "0 ";
            } else {
                for (int k = 0; k < n_cells; ++k) {
                    ofs << field[start_idx + k] << " ";
                }
            }
            ofs << "\n        </DataArray>\n";
        };

        write_scalar_field("Pressure",    surfphys_soa.pressure);
        write_scalar_field("Temperature", surfphys_soa.temperature);
        write_scalar_field("Tau1",        surfphys_soa.tau1);
        write_scalar_field("Tau2",        surfphys_soa.tau2);
        write_scalar_field("dTdn",        surfphys_soa.dTdn);

        write_int_field("Rank",       surfphys_soa.rank);
        write_int_field("Level",      surfphys_soa.lev);
        write_int_field("IP_quality", surfphys_soa.ip_quality);

        ofs << "      </CellData>\n";
        ofs << "    </Piece>\n";
        ofs << "  </PolyData>\n";
        ofs << "</VTKFile>\n";

        ofs.close();
    }
}

// ============================================================================
// 3. buildCSR
// ============================================================================

private:

AMREX_FORCE_INLINE
void buildCSR()
{
    const int myrank = amrex::ParallelDescriptor::MyProc();
    const int nlevel = amr_p->finestLevel() + 1;

    faces_per_level.clear();
    faces_per_level.resize(nlevel);

    Vector<int> nfab_per_level(nlevel);

    // ------------------------------------------------------------
    // 0) Initialization: Prepare CSR structure for each level
    // ------------------------------------------------------------
    for (int lev = 0; lev < nlevel; ++lev) {
        const int nfab_local = bmf_a[lev]->local_size();
        nfab_per_level[lev] = nfab_local;
        auto& csr = faces_per_level[lev];
        csr.fab_offsets.resize(nfab_local + 1, 0);
    }

    // ------------------------------------------------------------
    // 1) Pass-1: Count faces per FAB
    // ------------------------------------------------------------
    for (int f = 0; f < ntotalfaces; ++f) {
        if (surfphys_soa.elemfound[f] != 1) continue;
        if (surfphys_soa.rank[f] != myrank) continue;

        const int lev = surfphys_soa.lev[f];
        if (lev < 0 || lev >= nlevel) continue;

        const int ifab = surfphys_soa.ifab[f];
        if (ifab < 0 || ifab >= nfab_per_level[lev]) continue;

        faces_per_level[lev].fab_offsets[ifab]++;
    }

    // ------------------------------------------------------------
    // 2) Pass-2: Prefix Sum (Compute Offsets) & Allocation
    // ------------------------------------------------------------
    Vector<Vector<int>> cursor(nlevel);

    for (int lev = 0; lev < nlevel; ++lev) {
        auto& csr = faces_per_level[lev];
        const int nfab_local = nfab_per_level[lev];

        cursor[lev].resize(nfab_local);

        int offset = 0;
        for (int ifab = 0; ifab < nfab_local; ++ifab) {
            int cnt = csr.fab_offsets[ifab];
            csr.fab_offsets[ifab] = offset;
            cursor[lev][ifab]     = offset;
            offset += cnt;
        }
        csr.fab_offsets[nfab_local] = offset;
        csr.face_indices.resize(offset);
    }

    // ------------------------------------------------------------
    // 3) Pass-3: Fill face indices
    // ------------------------------------------------------------
    for (int f = 0; f < ntotalfaces; ++f) {
        if (surfphys_soa.elemfound[f] != 1) continue;
        if (surfphys_soa.rank[f] != myrank) continue;

        const int lev = surfphys_soa.lev[f];
        if (lev < 0 || lev >= nlevel) continue;

        const int ifab = surfphys_soa.ifab[f];
        if (ifab < 0 || ifab >= nfab_per_level[lev]) continue;

        const int pos = cursor[lev][ifab]++;
        faces_per_level[lev].face_indices[pos] = f;
    }

    // Update filled_elems count
    surfphys_soa.filled_elems = 0;
    for (int lev = 0; lev < nlevel; ++lev) {
        auto& csr = faces_per_level[lev];
        surfphys_soa.filled_elems += csr.face_indices.size();
    }
}

// ============================================================================
// 4. getGeomIdx
// ============================================================================

AMREX_FORCE_INLINE
int getGeomIdx(int elemIdx) const
{
    if (elemIdx < 0 || elemIdx >= ntotalfaces) return -1;

    for (int i = 0; i < ngeom; ++i) {
        if (elemIdx >= geom_offsets[i] && elemIdx < geom_offsets[i + 1]) {
            return i;
        }
    }
    return -1;
}

// ============================================================================
// 5. read_geom
// ============================================================================

void read_geom()
{
    ParmParse pp;
    Vector<std::string> files_a;
    bool plot_surf = false;

    pp.getarr("ib.filename", files_a);
    pp.query("ib.plot_surf", plot_surf);
    int skip_validation = 0;
    pp.query("ib.skip_validation", skip_validation);

    if (files_a.empty()) {
        amrex::Warning("ib.filename is empty: no geometry files provided");
    }

    // Release inside/outside testers
    for (auto* p : inout_fa) { delete p; }
    inout_fa.clear();

    // Resize all containers to match new geometry count
    this->ngeom = static_cast<int>(files_a.size());
    this->geom_a.resize(this->ngeom);
#ifdef AMREX_USE_CGAL
    this->tree_a.resize(this->ngeom);
    this->idxmap_a.resize(this->ngeom);
#else
    this->bvh_a.resize(this->ngeom);
#endif
    this->bbox_a.resize(this->ngeom);
    this->LocalFrame_a.clear();
    this->SurfElem_a.clear();
    this->geom_offsets.resize(this->ngeom + 1);
    this->inout_fa.resize(this->ngeom);
    this->ntotalfaces = 0;

    // Store geometry names (strip path and extension from filenames)
    this->geom_names.resize(this->ngeom);
    for (int i = 0; i < this->ngeom; ++i) {
        std::string name = files_a[i];
        // Strip directory path
        auto slash = name.find_last_of("/\\");
        if (slash != std::string::npos) name = name.substr(slash + 1);
        // Strip extension
        auto dot = name.find_last_of('.');
        if (dot != std::string::npos) name = name.substr(0, dot);
        this->geom_names[i] = name;
    }

    // Determine minimum dx across all levels for setting polygon tolerance
    Real min_dx = std::numeric_limits<Real>::max();
    for (int lev = 0; lev < dx_a.size(); ++lev) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            min_dx = std::min(min_dx, dx_a[lev][d]);
        }
    }
    min_dx *= Real(0.5);

    for (int i = 0; i < ngeom; i++) {
        Print() << "----------------------------------" << std::endl;

#if (AMREX_SPACEDIM == 2)

        if (!read_polygon_2d(files_a[i], geom_a[i], min_dx)) {
            amrex::Abort(std::string("Invalid 2D geometry filename: ") + files_a[i]);
        }
        Print() << "Geometry (i=" << i << ") " << files_a[i] << " read" << std::endl;
        Print() << "Number of vertices in polygon: " << geom_a[i].size() << "\n";

#ifdef AMREX_USE_CGAL
        tree_a[i].insert(geom_a[i].edges_begin(), geom_a[i].edges_end());
        tree_a[i].build();
        Print() << "CGAL AABB tree constructed" << std::endl;

        inout_fa[i] = new inside_t(geom_a[i]);
        Print() << "2D in/out testing functor constructed for polygon " << files_a[i] << "\n";
#else
        bvh_a[i].build(geom_a[i]);
        Print() << "BVH constructed" << std::endl;

        inout_fa[i] = new inside_t(geom_a[i], bvh_a[i]);
        Print() << "2D in/out testing functor constructed for polygon " << files_a[i] << "\n";
#endif

        bbox_a[i] = geom_a[i].bbox();
        ntotalfaces += geom_a[i].size();

#elif (AMREX_SPACEDIM == 3)

#ifdef AMREX_USE_CGAL
        if (!PMP::IO::read_polygon_mesh(files_a[i], geom_a[i])) {
            amrex::Abort(std::string("Invalid geometry filename: ") + files_a[i]);
        }
#else
        if (!read_mesh(files_a[i], geom_a[i])) {
            amrex::Abort(std::string("Invalid geometry filename: ") + files_a[i]);
        }
#endif
        Print() << "Geometry (i=" << i << ") " << files_a[i] << " read" << std::endl;
        Print() << "Number of facets " << geom_a[i].size_of_facets() << std::endl;

        if (!geom_a[i].is_closed()) {
            if (skip_validation) {
                amrex::Warning("IBM mesh is not closed (watertight) — validation skipped by ib.skip_validation=1");
            } else {
                amrex::Abort("IBM mesh validation failed: Mesh is not closed.");
            }
        }

#ifdef AMREX_USE_CGAL
        tree_a[i].insert(faces(geom_a[i]).first, faces(geom_a[i]).second, geom_a[i]);
        tree_a[i].build();
        Print() << "CGAL AABB tree constructed" << std::endl;

        inout_fa[i] = new inside_t(geom_a[i]);
        Print() << "In out testing function constructed for geometry " << files_a[i] << "\n";
#else
        bvh_a[i].build(geom_a[i]);
        Print() << "BVH constructed" << std::endl;

        inout_fa[i] = new inside_t(geom_a[i], bvh_a[i]);
        Print() << "In out testing function constructed for geometry " << files_a[i] << "\n";
#endif

#ifdef AMREX_USE_CGAL
        bbox_a[i] = PMP::bbox(geom_a[i]);
#else
        bbox_a[i] = geom_a[i].bbox();
#endif
        ntotalfaces += geom_a[i].size_of_facets();

#endif

        this->geom_offsets[i] = static_cast<int>(this->LocalFrame_a.size());
#ifdef AMREX_USE_CGAL
        build_geometry_cache(geom_a[i], SurfElem_a, LocalFrame_a,
                             idxmap_a[i], this->geom_offsets[i], i);
#else
        build_geometry_cache(geom_a[i], SurfElem_a, LocalFrame_a, this->geom_offsets[i], i);
#endif
    } // end loop over geometries

    this->geom_offsets[ngeom] = static_cast<int>(this->LocalFrame_a.size());

    if ((static_cast<int>(SurfElem_a.size()) != ntotalfaces) ||
        (static_cast<int>(LocalFrame_a.size()) != ntotalfaces)) {
        amrex::Print() << "Warning: Geometry cache size mismatch (degenerate faces skipped)\n"
                       << "  Expected total faces     : " << ntotalfaces << "\n"
                       << "  Actual SurfElem_a size   : " << SurfElem_a.size() << "\n"
                       << "  Actual LocalFrame_a size : " << LocalFrame_a.size() << "\n";
        ntotalfaces = static_cast<int>(SurfElem_a.size());
    }

    check_ibm_geometry_consistency(ngeom, geom_a.data(), inout_fa.data(), files_a.data());

    if (interior_is_solid) {
        Print() << "Interior of geometry is marked as SOLID" << std::endl;
    } else {
        convert_inout(LocalFrame_a);
        Print() << "Interior of geometry is marked as FLUID" << std::endl;
    }

    Print() << "----------------------------------" << std::endl;
    Print() << "----------------------------------" << std::endl;

    if (plot_surf) {
        Print() << "Total number of faces across all geometries: " << ntotalfaces << std::endl;
        surfimp_soa.resize(ntotalfaces);
        surfphys_soa.resize(ntotalfaces);
    }
}

#endif // EIB_IO_H_
