#!/bin/bash

for i in 16 32 64 128
do
    echo "\nRunning MMS test N = $i\n"
    mpirun -np 16 ./Cerisse3d.gnu.MPI.ex inputs_mms amr.n_cell=$i $i $i amr.data_log=euler/teno5-"$i".log cns.do_visc=0 prob.mach=1.0 cns.recon_scheme=6 cns.use_hybrid_scheme=0 cns.cfl=0.3 # cns.rk_order=2 cns.amr_interp_order=4
    # mpirun -np 16 ./Cerisse3d.gnu.MPI.ex inputs_mms amr.n_cell=$i $i $i amr.data_log=ns/teno5-"$i".log cns.do_visc=1 prob.reynolds=1.0 prob.prandtl=1.0 prob.mach=1.0 cns.recon_scheme=6 cns.use_hybrid_scheme=0 cns.cfl=0.3 # cns.rk_order=2 cns.amr_interp_order=4
    # mpirun -np 16 ./Cerisse3d.gnu.MPI.ex inputs_mms amr.n_cell=$i $i $i amr.data_log=amr4/teno5-"$i".log amr.max_level=1 cns.do_visc=1 prob.reynolds=1.0 prob.prandtl=1.0 prob.mach=1.0 cns.recon_scheme=6 cns.use_hybrid_scheme=0 cns.cfl=0.3 cns.amr_interp_order=4
done