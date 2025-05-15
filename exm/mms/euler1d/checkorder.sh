#!/bin/bash

# List of meshes to run
MESH=( 
	"16"
	"32"
	"64"
	"128"
	"256"
)

RUN="./main1d.gnu.ex"

echo "Checking  Spatial Order Euler ..."
echo "  "
echo "  "

# run and create directories 
for grid in "${MESH[@]}"; do
   $RUN input_dir/inputs$grid
    mv plot plot$grid
    printf " run  mesh with %s nodes \n "  "$grid"
done

# plot order (optional)
#python checkorder.py
