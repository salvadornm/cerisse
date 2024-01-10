#!/bin/bash

for i in 16 32 64 128
do
    echo "Running MMS test N = $i"
    mpiexec -np 16 ./Cerisse3d.gnu.TEST.MPI.ex inputs amr.n_cell=$i $i $i amr.data_log=wenoz3-"$i".log cns.recon_scheme=3
    mpiexec -np 16 ./Cerisse3d.gnu.TEST.MPI.ex inputs amr.n_cell=$i $i $i amr.data_log=wenoz5-"$i".log cns.recon_scheme=5
done