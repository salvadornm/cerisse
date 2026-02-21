
#!/bin/bash

# Colors
GREEN="\033[0;32m"
RED="\033[0;31m"
NC="\033[0m" # No color

# List of dirs to compile
MAKE_DIRS=( 
	"../exm/numerics/covo/"
	"../exm/numerics/riemann/"
	"../exm/numerics/shuosher/"
	"../exm/auto_ignition/"
	"../exm/reactive_sod/"
	"../exm/shock_reflect/"
	"../exm/ebm/forward_step/"
	"../exm/viscwall/"
	"../exm/ebm/cylinder/"
	"../exm/ebm/cylinder_visc/"
	"../exm/ebm/combustor/"
	"../exm/planar_flame/"	
	"../exm/ibm/sphere/"
	"../exm/ibm/srp/"
	"../exm/mms/navsto3d/"
	"../tst/tutorial/"
)
# List of dirs that include reaction
dirs_reacc=("../exm/auto_ignition/" "../exm/reactive_sod/" "../exm/planar_flame/" "../exm/ibm/srp/")
	
# make instruction
MAKECOMP="make -j8"

echo "Checking  Cerisse Examples (be patient) ..."
echo "  "
echo "  "

#header
printf "%-30s\t%-15s\t%-15s\n" "Directory" "Clean" "Build"
printf "%-30s\t%-15s\t%-15s\n" "---------" "-----" "-----"

for dir in "${MAKE_DIRS[@]}"; do
    
	clean_result=" "
	build_result=" "

	#echo -n "Cleaning in $dir..."
    	if make -C "$dir" clean &> /dev/null; then
		#echo -e "                      ${GREEN}[SUCCESS]${NC}"
		clean_result="${GREEN}[SUCCESS]${NC}"
	else 
		#echo -e "                      ${RED}[FAILED]${NC}   "
		clean_result="${RED}[FAILED]${NC}"
	fi

    # If dir is in dirs_reacc, build TPL first
	if [[ " ${dirs_reacc[@]} " =~ " ${dir} " ]]; then
		make -C "$dir" TPLclean &> /dev/null;
		echo  "Building Chemistry in $dir..."
		if ! make -C "$dir" TPL > build.log 2>&1; then
			echo "===== TPL build failed in $dir ====="
    		tail -n 10 build.log
    		exit 1
			echo "===================================="
		fi		    	
	fi 

	#echo -n "Running make in $dir..."
	if $MAKECOMP -C "$dir" > build.log 2>&1; then
		#echo -e "                      ${GREEN}[PASS]${NC}"
		build_result="${GREEN}[PASS]${NC}" 
	else 
		#echo -e "                      ${RED}[FAIL]${NC}"
		build_result="${RED}[FAIL]${NC}"
		cat build.log
		echo "===== compile failed in $dir ====="
    	tail -n 10 build.log
    	exit 1
		echo "===================================="
	fi
	# remove compiled files, executables and log files
	rm -rf "${dir:?}/tmp_build_dir"
	rm -rf "${dir:?}/*.ex"
	rm  *.log

	printf "%-30s\t%-15b\t%-15b\n" "$dir" "$clean_result" "$build_result"
done		
printf "%-30s\t%-15s\t%-15s\n" "---------" "-----" "-----"
