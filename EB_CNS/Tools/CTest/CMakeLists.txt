#============================================================================
# BUILDING TEST CASES
# 
# Steps for addding a new test named new_test from the example NEW_CASE:
# 1. Create an input file in the directory NEW_CASE, named new_test.ini.
# 2. Add the following line in the end of this file: 
#     add_test_r(new_test NEW_CASE labels)
#    where "labels" specifies the group of the test case.
#
# TODO: Enable comparison with benchmark (Gold) pltfiles using fcompare.
#============================================================================

set(CTEST_CUSTOM_MAXIMUM_FAILED_TEST_OUTPUT_SIZE "262144")
set(CTEST_CUSTOM_MAXIMUM_PASSED_TEST_OUTPUT_SIZE "51200")

# if(CNS_ENABLE_FCOMPARE_FOR_TESTS)
#   if("${CNS_REFERENCE_GOLDS_DIRECTORY}" STREQUAL "")
#     message(FATAL_ERROR "To reference gold files, CNS_REFERENCE_GOLDS_DIRECTORY must be set and exist")
#   else()
#     if(EXISTS ${CNS_REFERENCE_GOLDS_DIRECTORY})
#       set(GOLD_FILES_DIRECTORY ${CNS_REFERENCE_GOLDS_DIRECTORY}/${CMAKE_SYSTEM_NAME}/${CMAKE_CXX_COMPILER_ID}/${CMAKE_CXX_COMPILER_VERSION})
#       message(STATUS "Test golds directory for fcompare: ${GOLD_FILES_DIRECTORY}")
#     else()
#       message(FATAL_ERROR "Specified directory for reference gold files does not exist: ${CNS_REFERENCE_GOLDS_DIRECTORY}")
#     endif()
#   endif()
# endif()

# if(CNS_SAVE_GOLDS)
#   if("${CNS_SAVED_GOLDS_DIRECTORY}" STREQUAL "")
#     message(FATAL_ERROR "To save gold files, CNS_SAVED_GOLDS_DIRECTORY must be set and exist")
#   else()
#     if(EXISTS ${CNS_SAVED_GOLDS_DIRECTORY})
#       set(SAVED_GOLDS_DIRECTORY ${CNS_SAVED_GOLDS_DIRECTORY}/${CMAKE_SYSTEM_NAME}/${CMAKE_CXX_COMPILER_ID}/${CMAKE_CXX_COMPILER_VERSION})
#       message(STATUS "Gold files will be saved to: ${SAVED_GOLDS_DIRECTORY}")
#     else()
#       message(FATAL_ERROR "Specified directory for saving gold files does not exist: ${CNS_SAVED_GOLDS_DIRECTORY}")
#     endif()
#   endif()
# endif()

######################### Functions for adding tests #########################
macro(setup_test)
  # Set variables for respective binary and source directories for the test
  set(CURRENT_TEST_SOURCE_DIR ${CMAKE_SOURCE_DIR}/EB_CNS/Exec/${TEST_EXE_DIR})
  set(CURRENT_TEST_BINARY_DIR ${CMAKE_BINARY_DIR}/EB_CNS/Exec/${TEST_EXE_DIR}/tests/${TEST_NAME})
  set(CURRENT_TEST_EXE ${CMAKE_BINARY_DIR}/EB_CNS/Exec/${TEST_EXE_DIR}/${TEST_EXE_DIR}.ex)
  # Gold files should be submodule organized by machine and compiler (these are output during configure)
  # set(PLOT_GOLD ${GOLD_FILES_DIRECTORY}/${TEST_EXE_DIR}/tests/${TEST_NAME}/plt00010)
  # Test plot is currently expected to be after 10 steps
  # set(PLOT_TEST ${CURRENT_TEST_BINARY_DIR}/plt00010)
  # Find fcompare
  # if(CNS_ENABLE_FCOMPARE_FOR_TESTS)
  #   set(FCOMPARE ${CMAKE_BINARY_DIR}/Submodules/AMReX/Tools/Plotfile/fcompare)
  # endif()
  # Make working directory for test
  file(MAKE_DIRECTORY ${CURRENT_TEST_BINARY_DIR})
  # Gather all files in source directory for test
  # file(GLOB TEST_FILES "${CURRENT_TEST_SOURCE_DIR}/*.dat" "${CURRENT_TEST_SOURCE_DIR}/*.py")
  # Copy files to test working directory
  file(COPY ${CURRENT_TEST_SOURCE_DIR}/${TEST_NAME}.ini DESTINATION "${CURRENT_TEST_BINARY_DIR}/")
  # file(COPY ${TEST_FILES} DESTINATION "${CURRENT_TEST_BINARY_DIR}/")
  # Set some default runtime options for all tests
  set(RUNTIME_OPTIONS "amr.plot_file=plt amr.checkpoint_files_output=0 amr.plot_files_output=1") #amrex.the_arena_is_managed=0
  # Trap floating point errors (fpe)
  set(RUNTIME_OPTIONS "${RUNTIME_OPTIONS} amrex.signal_handling=1 amrex.fpe_trap_invalid=1 amrex.fpe_trap_zero=1 amrex.fpe_trap_overflow=1")

  if(CNS_ENABLE_MPI)
    if(NOT DEFINED CNS_NP)
      set(CNS_NP 4)
    endif()
    set(MPI_COMMANDS "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${CNS_NP} ${MPIEXEC_PREFLAGS}")
  else()
    set(CNS_NP 1)
    unset(MPI_COMMANDS)
  endif()
  # Use fcompare to test diffs in plots against gold files
  # if(CNS_ENABLE_FCOMPARE_FOR_TESTS AND (NOT "${TEST_NAME}" MATCHES "hdf5$"))
  #   if(CNS_ENABLE_CUDA)
  #     set(FCOMPARE_TOLERANCE "-r 1e-12 --abs_tol 1.0e-12")
  #   endif()
  #   set(FCOMPARE_COMMAND "&& ${MPI_COMMANDS} ${FCOMPARE} ${FCOMPARE_TOLERANCE} ${PLOT_TEST} ${PLOT_GOLD}")
  # endif()
  # if(CNS_SAVE_GOLDS)
  #   file(MAKE_DIRECTORY ${SAVED_GOLDS_DIRECTORY}/${TEST_EXE_DIR}/${TEST_NAME})
  #   set(SAVE_GOLDS_COMMAND "&& cp -R ${PLOT_TEST} ${SAVED_GOLDS_DIRECTORY}/${TEST_EXE_DIR}/${TEST_NAME}/")
  # endif()
endmacro(setup_test)

# Standard regression test
function(add_test_r TEST_NAME TEST_EXE_DIR LABEL)
  setup_test()
  # set(RUNTIME_OPTIONS "max_step=10 ${RUNTIME_OPTIONS}")
  add_test(${TEST_NAME} 
           sh -c "${MPI_COMMANDS} ${CURRENT_TEST_EXE} ${MPIEXEC_POSTFLAGS} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.ini ${RUNTIME_OPTIONS} ${EXTRA_RUNTIME_OPTIONS} > ${TEST_NAME}.log ${SAVE_GOLDS_COMMAND} ${FCOMPARE_COMMAND}")
  set_tests_properties(${TEST_NAME} 
                       PROPERTIES 
                        TIMEOUT 18000 
                        PROCESSORS ${CNS_NP} 
                        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/" 
                        LABELS "${LABEL}"
                        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log")
endfunction(add_test_r)

############################# Regression tests #############################

add_test_r(moving-vortex ConvectingVortex  "Reg;2D")
add_test_r(flat-plate    FlatPlate         "Reg;2D;Diffuse")
add_test_r(mixlayer      MixingLayer       "Reg;2D;React;Diffuse")
add_test_r(react-sod     ReactiveShockTube "Reg;2D;React")
add_test_r(dmr           ShockReflect      "Reg;2D")
add_test_r(supersonic-jet SupersonicJet    "Reg;3D;React")
add_test_r(tgv-64        TGV               "Reg;3D;Diffuse")

if(CNS_ENABLE_EB)
  add_test_r(ddt               DDT             "Reg;2D;React;EB;Diffuse")
  add_test_r(detonation-diffrac DetonDiffrac   "Reg;2D;React;EB")
  add_test_r(flat-plate-eb     FlatPlate_EB    "Reg;2D;EB;Diffuse")
  add_test_r(forward-step      ForwardStep     "Reg;2D;EB")
  add_test_r(dmr-eb            ShockReflect_EB "Reg;2D;EB")
  add_test_r(double-ramp-isothermal DoubleRamp "Reg;3D;EB;Diffuse")
endif()

############################ Verification tests ############################

# ENABLE_MASA
# add_test_v1(symmetry MMS)
# add_test_v1(eb-symmetry EB-C1)

################################ Unit tests ################################

#add_test_u(unit-tests)