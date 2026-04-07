#!/bin/bash

# Usage: ./pp_modify.sh

# cd "$(dirname "$0")" || exit 1

# for PELEPHYSICS25

echo "... Modifying PelePhysics/Source/PelePhysicsConstraints.H"
FILE="./lib/PelePhysics/Source/PelePhysicsConstraints.H"
sed -i '39,43{
  s|^[[:space:]]*//|//|;       # Keep already-commented lines as is
  t done
  s|^|//|
  :done
}' "$FILE"

echo "... Modifying PelePhysics/Source/Reactions/ReactorCvodePreconditioner.cpp"
FILE="./lib/PelePhysics/Source/Reactions/ReactorCvodePreconditioner.cpp"
sed -i '115{/cuS_st/{ # Keep already-commented lines as is
  s/.*/#ifdef AMREX_USE_FLOAT\
    cuS_st = cusolverSpScsrqrsvBatched\
#else\
    cuS_st = cusolverSpDcsrqrsvBatched\
#endif\
   (/
}}' "$FILE"

echo "... Modifying PelePhysics/Source/Reactions/ReactorUtils.cpp"
FILE="./lib/PelePhysics/Source/Reactions/ReactorUtils.cpp"
sed -i '31{
  s|^[[:space:]]*//|//|;       # Keep already-commented lines as is
  t done
  s|^|//|
  :done
}' "$FILE"

echo "... Modifying PelePhysics/Source/Reactions/ReactorCvode.cpp"
FILE="./lib/PelePhysics/Source/Reactions/ReactorCvode.cpp"
sed -i '1164{/cusolver_status/{ # Keep already-commented lines as is
  s/.*/#ifdef AMREX_USE_FLOAT\
    cusolver_status = cusolverSpScsrqrBufferInfoBatched\
#else\
    cusolver_status = cusolverSpDcsrqrBufferInfoBatched\
#endif\
    (/
}}' "$FILE"

echo "... Modifying PelePhysics/Source/Reactions/ReactorCvode.cpp (Adding cleanup)"
FILE="./lib/PelePhysics/Source/Reactions/ReactorCvode.cpp"
sed -i '1355{/return (1);/{
  s|.*|    //clean up even when error! \
    N_VDestroy(y);\
    CVodeFree(\&cvode_mem);\
    if (LS != nullptr) {\
      SUNLinSolFree(LS);\
    }\
    if (NLS != nullptr) {\
      SUNNonlinSolFree(NLS);\
    }\
    if (A != nullptr) {\
      SUNMatDestroy(A);\
    }\
    freeUserData(udata);\
    return (1);|
}}' "$FILE"

echo "DONE"