#!/bin/bash

# Convert addresses in AMReX Backtrace file to human-readable lines
# Usage: ./convert-backtrace.sh [EXECUTABLE] [BACKTRACE_FILE]

file_name=$2

while read line; do

  address=$(echo $line | grep -Eo "\(\+.*?\)" | sed 's/^..//g; s/.$//g')
  addrlen=$(echo $address | wc -c)

  # output=$(addr2line -Cpfie $1 $address)
  if [ $addrlen -gt 1 ]; then
    echo "$(addr2line -Cpfie $1 $address)"
  fi

done < $file_name