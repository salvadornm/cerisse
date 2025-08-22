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

echo "DONE"