#!/bin/bash

# Auto-format all *.cpp and *.H files using clang-format.
# See the .clang-format file for detailed formatting rules.
# Usage: run this from the cerisse/EB_CNS directory, i.e. 
#        cd ../.. && ./Tools/clang-format/format-all.sh

# Find all files in all subdirectories
files=$(find Source Exec Tools -type f -name "*.cpp" -or -name "*.H")

# For each file, format it
for file in $files; do
  echo "Formatting $file..."
  clang-format -i $file
done