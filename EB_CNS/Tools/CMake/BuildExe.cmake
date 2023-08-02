#============================================================================
# BUILDING EXECUTABLE FOR EXAMPLE CASES
# 
# Follow these steps for addding a new case named NEW_CASE:
# 1. Put all problem files (prob_parm.H, prob.H, prob.cpp) in Exec/NEW_CASE.
# 2. Create a CMakeLists.txt in NEW_CASE, which includes all compile options
#    and EOS, Mechanism, and Transport options. At the end of the file, call 
#    include(BuildExe) to call this file.
# 3. In Exec/CMakeLists.txt, add the line add_subdirectory(NEW_CASE).
#============================================================================

get_filename_component(DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(exe_name ${DIR_NAME}.ex) #-${CNS_EOS_MODEL}-${CNS_CHEMISTRY_MODEL}-${CNS_TRANSPORT_MODEL}

add_executable(${exe_name} "")

target_sources(${exe_name}
  PRIVATE
    prob_parm.H
    prob.H
    prob.cpp
)

############################# Include PelePhysics #############################
target_include_directories(${exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

set(PELE_PHYSICS_SRC_DIR ${CMAKE_SOURCE_DIR}/Submodules/PelePhysics)
set(PELE_PHYSICS_BIN_DIR ${CMAKE_BINARY_DIR}/Submodules/PelePhysics/${exe_name})
target_include_directories(${exe_name} SYSTEM PRIVATE "${PELE_PHYSICS_SRC_DIR}/Source")

# Transport
set(PP_TRANSPORT_DIR "${PELE_PHYSICS_SRC_DIR}/Transport")
target_include_directories(${exe_name} SYSTEM PRIVATE ${PP_TRANSPORT_DIR})
target_sources(${exe_name} PRIVATE
    ${PP_TRANSPORT_DIR}/Transport.H
    ${PP_TRANSPORT_DIR}/Transport.cpp
    ${PP_TRANSPORT_DIR}/TransportParams.H
    ${PP_TRANSPORT_DIR}/TransportTypes.H
    ${PP_TRANSPORT_DIR}/Simple.H
    ${PP_TRANSPORT_DIR}/Sutherland.H
    ${PP_TRANSPORT_DIR}/Constant.H
  )
if (CNS_TRANSPORT_MODEL STREQUAL "Simple")
  target_compile_definitions(${exe_name} PRIVATE USE_SIMPLE_TRANSPORT)
elseif(CNS_TRANSPORT_MODEL STREQUAL "Constant")
  target_compile_definitions(${exe_name} PRIVATE USE_CONSTANT_TRANSPORT)
elseif(CNS_TRANSPORT_MODEL STREQUAL "Sutherland")
  target_compile_definitions(${exe_name} PRIVATE USE_SUTHERLAND_TRANSPORT)
endif()

# EOS
set(PP_EOS_DIR "${PELE_PHYSICS_SRC_DIR}/Eos")
target_include_directories(${exe_name} SYSTEM PRIVATE ${PP_EOS_DIR})
target_sources(${exe_name} PRIVATE
               ${PP_EOS_DIR}/EOS.cpp
               ${PP_EOS_DIR}/EOS.H
               ${PP_EOS_DIR}/Fuego.H
               ${PP_EOS_DIR}/GammaLaw.H
               ${PP_EOS_DIR}/SRK.H)
if (CNS_EOS_MODEL STREQUAL "Fuego")
  target_compile_definitions(${exe_name} PRIVATE USE_FUEGO_EOS)
elseif(CNS_EOS_MODEL STREQUAL "GammaLaw")
  target_compile_definitions(${exe_name} PRIVATE USE_GAMMALAW_EOS)
elseif(CNS_EOS_MODEL STREQUAL "Soave-Redlich-Kwong")
  target_compile_definitions(${exe_name} PRIVATE USE_SRK_EOS)
endif()

# Mechanism
if(NOT DEFINED PP_MECHANISM_DIR)
  set(PP_MECHANISM_DIR "${PELE_PHYSICS_SRC_DIR}/Support/Mechanism/Models/${CNS_CHEMISTRY_MODEL}")
endif()
target_include_directories(${exe_name} SYSTEM PRIVATE ${PP_MECHANISM_DIR})
target_include_directories(${exe_name} SYSTEM PRIVATE ${PELE_PHYSICS_SRC_DIR}/Support/Mechanism/Models)
target_sources(${exe_name} PRIVATE ${PP_MECHANISM_DIR}/mechanism.cpp 
                                   ${PP_MECHANISM_DIR}/mechanism.H)

# Reactor
target_include_directories(${exe_name} PRIVATE ${PELE_PHYSICS_SRC_DIR}/Reactions)
target_sources(${exe_name}
  PRIVATE
    ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorArkode.H
    ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorArkode.cpp
    ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorBase.H
    ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorBase.cpp
    ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorCvode.H
    ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorCvode.cpp
    ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorCvodeCustomLinSolver.H
    ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorCvodeCustomLinSolver.cpp
    ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorCvodeJacobian.H
    ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorCvodeJacobian.cpp
    ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorCvodePreconditioner.H
    ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorCvodePreconditioner.cpp
    ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorCvodeUtils.H
    ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorCvodeUtils.cpp
    ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorNull.H
    ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorNull.cpp
    ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorRK64.H
    ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorRK64.cpp
    ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorTypes.H
    ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorUtils.H
    ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorUtils.cpp
)

# Sundials
target_include_directories(${exe_name} PRIVATE ${AMREX_DIR}/Src/Extern/SUNDIALS)
target_sources(${exe_name} PRIVATE ${AMREX_DIR}/Src/Extern/SUNDIALS/AMReX_Sundials.H
                                    ${AMREX_DIR}/Src/Extern/SUNDIALS/AMReX_Sundials_Core.cpp
                                    ${AMREX_DIR}/Src/Extern/SUNDIALS/AMReX_Sundials_Core.H
                                    ${AMREX_DIR}/Src/Extern/SUNDIALS/AMReX_NVector_MultiFab.cpp
                                    ${AMREX_DIR}/Src/Extern/SUNDIALS/AMReX_NVector_MultiFab.H
                                    ${AMREX_DIR}/Src/Extern/SUNDIALS/AMReX_SUNMemory.cpp
                                    ${AMREX_DIR}/Src/Extern/SUNDIALS/AMReX_SUNMemory.H)
target_link_libraries(${exe_name} PRIVATE sundials_arkode sundials_cvode)

if(CNS_ENABLE_CUDA)
  target_link_libraries(${exe_name} PRIVATE sundials_nveccuda sundials_sunlinsolcusolversp sundials_sunmatrixcusparse)
elseif(CNS_ENABLE_HIP)
  target_link_libraries(${exe_name} PRIVATE sundials_nvechip)
elseif(CNS_ENABLE_SYCL)
  target_link_libraries(${exe_name} PRIVATE sundials_nvecsycl)
endif()

############################# Include CNS Source #############################

set(SRC_DIR ${CMAKE_SOURCE_DIR}/EB_CNS/Source)
target_include_directories(${exe_name} PRIVATE ${SRC_DIR})
target_include_directories(${exe_name} PRIVATE ${SRC_DIR}/diffusion)
target_include_directories(${exe_name} PRIVATE ${SRC_DIR}/hydro)
target_include_directories(${exe_name} PRIVATE ${SRC_DIR}/pdf)
target_include_directories(${exe_name} PRIVATE ${SRC_DIR}/turb)
target_include_directories(${exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(${exe_name} PRIVATE ${CMAKE_BINARY_DIR})

target_sources(${exe_name}
    PRIVATE
      ${SRC_DIR}/CNS.H
      ${SRC_DIR}/derive.H
      ${SRC_DIR}/tagging.H
      ${SRC_DIR}/CNS_K.H
      ${SRC_DIR}/index_macros.H 
      ${SRC_DIR}/parm.H
      ${SRC_DIR}/default_parm.H
      ${SRC_DIR}/bcfill.H
      ${SRC_DIR}/bc_util.H
      ${SRC_DIR}/main.cpp
      ${SRC_DIR}/advance.cpp
      ${SRC_DIR}/advance_box.cpp
      ${SRC_DIR}/bcfill.cpp
      ${SRC_DIR}/derive.cpp
      ${SRC_DIR}/CNS.cpp
      ${SRC_DIR}/CNSBld.cpp
      ${SRC_DIR}/io.cpp
      ${SRC_DIR}/setup.cpp
      ${SRC_DIR}/react.cpp
      ${SRC_DIR}/diffusion/diffusion.H
      ${SRC_DIR}/hydro/hyperbolics.H 
      ${SRC_DIR}/hydro/recon.H
      ${SRC_DIR}/pdf/pdf_model.H 
      ${SRC_DIR}/pdf/pdf_model.cpp
      ${SRC_DIR}/pdf/random.H        
      ${SRC_DIR}/turb/LES.H 
      ${SRC_DIR}/turb/Smagorinsky.cpp
      ${SRC_DIR}/turb/WALE.cpp
      ${SRC_DIR}/turb/PaSR.H
)

if(CNS_ENABLE_EB)
  target_compile_definitions(${exe_name} PRIVATE CNS_USE_EB)

  target_sources(${exe_name}
    PRIVATE
    ${SRC_DIR}/custom_geometry.H
    ${SRC_DIR}/init_eb2.cpp
    ${SRC_DIR}/advance_box_eb.cpp
    ${SRC_DIR}/custom_geometry.cpp
    ${SRC_DIR}/diffusion/diffusion_eb.H
    ${SRC_DIR}/hydro/recon_eb.H
    ${SRC_DIR}/hydro/divop_eb.H
  )
endif()

############################# Link AMReX library #############################

target_link_libraries(${exe_name} PRIVATE AMReX::amrex)

if(CNS_ENABLE_CUDA)
  setup_target_for_cuda_compilation(${exe_name}) # this must be at the end
endif()