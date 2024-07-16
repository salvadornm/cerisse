#!/usr/bin/env bash

PELE_PHYSICS_HOME="../../../../Submodules/PelePhysics"
MECH_HOME="$(pwd)"

# 1. Convert CHEMKIN to YAML
# cd ${PELE_PHYSICS_HOME}/Support/ceptr
# poetry run ck2yaml --input=${MECH_HOME}/mechanism.inp --permissive

# 2. Convert YAML to mechanism
bash ${PELE_PHYSICS_HOME}/Support/Mechanism/Models/converter.sh -f "${MECH_HOME}/mechanism.yaml"
