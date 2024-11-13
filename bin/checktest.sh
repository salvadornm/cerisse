#!/bin/bash

if [ -n "$CERISSE_HOME" ]; then
    echo "CERISSE_HOME defined "
else
    echo "CERISSE_HOME undefined or empty"
    CERISSE_HOME=../..
fi

path_test1="$CERISSE_HOME/tst/test1"
path_ref="$CERISSE_HOME/bin/auxdat"

echo "-- Checking if code compiles ...  (be patient 3/4 mins)"

cd $path_test1

#clean & compile
make clean  &> /dev/null
make  -j4 &> /dev/null
compile_status=$?
if [ $compile_status -eq 0 ]; then
    echo "                              ...  [OK]"    
else
	echo "                 ... compilation [FAILED]"
	echo " ERROR check compilation output :"
    make
	exit 1
fi

# run

EXE=./main1d.gnu.MPI.ex

# run correctly

echo "-- Checking if  code runs ...  "
rm -rf time.log
$EXE inputs &> run.log


if [ $? -eq 0 ]; then
    echo "                              ...  [OK]"    
else
    echo "                           ... [FAILED]"
	echo " ERROR check run.log "
	exit 1
fi

#files (store reference in bin/auxfiles)
ref_file=$path_ref/timeref.dat
out_file=time.log

# results correct
echo "-- Checking if  code output is correct ...  "

diff --brief $out_file $ref_file >/dev/null
comp_value=$?

if [ $comp_value -eq 0 ]; then
    echo "                              ...  [OK]" 
    echo " Program ran successfully and output is correct."
    rm *.log
    rm -rf plot.* 
    rm -rf *.ex 
    rm -rf tmp_build_dir 
    exit 0  # Success
else    
    echo "                          ... [FAILED]"
    echo " got : "
    tail -n 1 $out_file
    echo " expected : "
    tail -n 1 $ref_file
    echo " ERROR check run.log "
	exit 1
fi


