#!/bin/bash

# This script installs MASA on the local machine

THIS_DIR="$(dirname "$(realpath "$0")")"

wget https://github.com/manufactured-solutions/MASA/releases/download/0.50.0/masa-0.50.0.tar.gz
tar xvfz masa-0.50.0.tar.gz
cd masa-0.50.0
./configure --"$THIS_DIR"/masa-0.50.0/install/ # this is the MASA install path, MASA_HOME in GNUmakefile
make
make check
make install