add_subdirectory(ConvectingVortex)
add_subdirectory(FlatPlate)
add_subdirectory(TGV)
add_subdirectory(MixingLayer)
add_subdirectory(ReactiveShockTube)
# add_subdirectory(Riemann2D) # which case?
add_subdirectory(ShockReflect)
add_subdirectory(SupersonicJet)

if(CNS_ENABLE_EB)
  add_subdirectory(DDT)
  add_subdirectory(DetonDiffrac)
  add_subdirectory(FlatPlate_EB)
  add_subdirectory(ForwardStep)
  add_subdirectory(ShockReflect_EB)
  add_subdirectory(DoubleRamp)
endif()

# Only available for 2D
# if(CNS_DIM STREQUAL 2)
#   ...
# endif()

# Only available for 3D
# if(CNS_DIM STREQUAL 3)
#   ...
# endif()