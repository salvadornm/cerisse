
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
	"../tst/tutorial/"
)
	
# make instruction
MAKECOMP="make -j4"

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
	#echo -n "Running make in $dir..."
	if $MAKECOMP -C "$dir" &> /dev/null; then
		#echo -e "                      ${GREEN}[PASS]${NC}"
		build_result="${GREEN}[PASS]${NC}"
	else 
		#echo -e "                      ${RED}[FAIL]${NC}"
		build_result="${RED}[FAIL]${NC}"
	fi

	printf "%-30s\t%-15b\t%-15b\n" "$dir" "$clean_result" "$build_result"
done		
printf "%-30s\t%-15s\t%-15s\n" "---------" "-----" "-----"
