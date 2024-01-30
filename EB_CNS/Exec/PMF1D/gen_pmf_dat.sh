#!/bin/bash
# Generate PMF data file using PelePhysics PMF utility
# It requires Python and poetry: install poetry using pip
# Usage: (OPTIONAL) PELE_PHYSICS_HOME=${PELE_PHYSICS_HOME} (OPTIONAL) OUTPUT_DIR=${OUTPUT_DIR} ./gen_pmf_dat.sh

if [ -z "$PELE_PHYSICS_HOME" ]; then
  PELE_PHYSICS_HOME=$(realpath "$(dirname "${BASH_SOURCE[0]}")"/../../Submodules/PelePhysics)
fi
if [ -z "$OUTPUT_DIR" ]; then
  OUTPUT_DIR=$(realpath .)
fi
echo "PELE_PHYSICS_HOME = $PELE_PHYSICS_HOME"
echo "OUTPUT_DIR        = $OUTPUT_DIR"

poetry -C $PELE_PHYSICS_HOME/Support/ceptr/ update
poetry -C $PELE_PHYSICS_HOME/Support/ceptr/ run python $PELE_PHYSICS_HOME/Utility/PMF/cantera_pmf_generator.py -pp $PELE_PHYSICS_HOME -m drm19 -f CH4 -d 0.016 -o $OUTPUT_DIR