#!/bin/bash

# Usage: ./pp_modify.sh

# cd "$(dirname "$0")" || exit 1

echo "... Modifying PelePhysics/Source/PelePhysicsConstraints.H"
FILE="../../../Submodules/PelePhysics/Source/PelePhysicsConstraints.H"
sed -i '30,35{
  s|^[[:space:]]*//|//|;       # Keep already-commented lines as is
  t done
  s|^|//|
  :done
}' "$FILE"

echo "... Modifying PelePhysics/Reactions/ReactorCvodeUtils.cpp"
FILE="../../../Submodules/PelePhysics/Reactions/ReactorCvodeUtils.cpp"
sed -i '16{
  s|^[[:space:]]*//|//|;       # Keep already-commented lines as is
  t done
  s|^|//|
  :done
}' "$FILE"

echo "... Modifying PelePhysics/Reactions/ReactorCvode.cpp"
FILE="../../../Submodules/PelePhysics/Reactions/ReactorCvode.cpp"
sed -i '1168{/cusolver_status/{ # Keep already-commented lines as is
  s/.*/#ifdef AMREX_USE_FLOAT\
    cusolver_status = cusolverSpScsrqrBufferInfoBatched\
#else\
    cusolver_status = cusolverSpDcsrqrBufferInfoBatched\
#endif\
    (/
}}' "$FILE"

echo "... Modifying PelePhysics/Reactions/ReactorCvodePreconditioner.cpp"
FILE="../../../Submodules/PelePhysics/Reactions/ReactorCvodePreconditioner.cpp"
sed -i '115{/cuS_st/{ # Keep already-commented lines as is
  s/.*/#ifdef AMREX_USE_FLOAT\
    cuS_st = cusolverSpScsrqrsvBatched\
#else\
    cuS_st = cusolverSpDcsrqrsvBatched\
#endif\
   (/
}}' "$FILE"

echo "DONE"