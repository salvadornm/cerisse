#!/usr/bin/env bash

MECH_HOME="$(pwd)"
MECH_FILE="${MECH_HOME}/mechanism.yaml"
bash ${PELE_PHYSICS_HOME}/Support/Mechanism/Models/converter.sh -f "${MECH_FILE}"
