#!/usr/bin/env bash

# Example usage: cd OConaire-9-21; ../convert.sh

PELE_PHYSICS_HOME="../../../../Submodules/PelePhysics/" # redefine this if needed
MECH_HOME="$(pwd)"

# Convert CHEMKIN files to mechanism.yaml 
cd ${PELE_PHYSICS_HOME}/Support/ceptr
poetry run ck2yaml --input=${MECH_HOME}/mechanism.inp --thermo=${MECH_HOME}/therm.dat --transport=${MECH_HOME}/tran.dat --permissive

# Generating mechanism.H and mechanism.cpp
cd ${MECH_HOME}
bash ${PELE_PHYSICS_HOME}/Support/Mechanism/Models/converter.sh -f "${MECH_HOME}"/mechanism.yaml

# Generate Make.package
echo "CEXE_headers+=mechanism.H
CEXE_sources+=mechanism.cpp" > Make.package