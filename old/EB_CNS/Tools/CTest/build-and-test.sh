#!/bin/bash

# Build and test solver using CMake + CTest
# Usage: run this from the root /cerisse directory, i.e.
#        cd ../../.. && ./EB_CNS/Tools/CTest/build-and-test.sh [2D/3D]

mkdir Build
cd Build

if [ "$1" == "2D" ]; then
  # 2D
  cmake -DCNS_DIM=2 \
        -DCNS_PRECISION=DOUBLE \
        -DCNS_ENABLE_TINY_PROFILE=ON \
        -DCNS_ENABLE_TESTING=ON \
        -DCNS_ENABLE_EB=ON \
        -DCNS_ENABLE_MPI=ON \
        -DCNS_NP=16 \
        -DEXTRA_RUNTIME_OPTIONS="max_step=100" \
        ..
elif [ "$1" == "3D" ]; then
  # 3D
  cmake -DCNS_DIM=3 \
        -DCNS_PRECISION=DOUBLE \
        -DCNS_ENABLE_TINY_PROFILE=ON \
        -DCNS_ENABLE_TESTING=ON \
        -DCNS_ENABLE_EB=ON \
        -DCNS_ENABLE_MPI=ON \
        -DCNS_NP=16 \
        -DEXTRA_RUNTIME_OPTIONS="max_step=10" \
        ..
else
  echo "Dimension of test must be 2D or 3D. Exit."  
  exit 1
fi

make -j16

ctest