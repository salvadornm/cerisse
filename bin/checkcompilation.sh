
#!/bin/bash

echo "Usage:"
echo "  ./checkcompilation.sh [MODE] [PELEPHYSICS=ON/OFF]"
echo ""
echo "Examples:"
echo "  ./checkcompilation.sh ALL"
echo "  ./checkcompilation.sh EBM"
echo "  ./checkcompilation.sh ALL OFF"
echo "  ./checkcompilation.sh STOP_IF_FAIL=TRUE"
echo ""

# Colors
GREEN="\033[0;32m"
RED="\033[0;31m"
NC="\033[0m" # No color

# List of dirs to compile
MAKE_DIRS=( 	
	"../tst/regtest/num/riemann/"
	"../tst/regtest/num/shuosher/"
	"../tst/regtest/num/covo/"
	"../tst/regtest/num/shock_reflect/"
	"../tst/regtest/visc/viscwall/"	
	"../tst/regtest/ebm/cylinder/"
	"../tst/regtest/ebm/forward_step/"
	"../tst/regtest/ebm/combustor_chem/"
	"../tst/regtest/ibm/sphere_cgal/"
	"../tst/regtest/ibm/polygon_bvh/"
	"../tst/regtest/ibm/srp_cgal_chem/"
	"../tst/regtest/react/autoigni_chem/"
	"../tst/regtest/react/sod_chem/"
	"../tst/regtest/react/flame_chem/"	
	"../tst/regtest/mms/navsto3d/"
	"../tst/regtest/tutorial/"
	"../tst/moving/airfoil/"
)

# List of dirs that include reaction
dirs_reacc=(
	"../tst/regtest/react/autoigni_chem/"
	"../tst/regtest/react/sod_chem/"
	"../tst/regtest/react/flame_chem/"	
	"../tst/regtest/ebm/combustor_chem/"
	"../tst/regtest/ibm/srp_cgal_chem/"
)

# -------------------------------
# Argument parsing
# -------------------------------
MODE=${1:-ALL}   		# default = ALL
PELEPHYSICS=${2:-ON}	# default = ON

FILTERED_DIRS=()

for dir in "${MAKE_DIRS[@]}"; do

    include_dir=true

    # ---- Mode filters ----
    if [[ "$MODE" == "EBM" ]]; then
        [[ "$dir" != *"/ebm/"* ]] && include_dir=false
    fi

    if [[ "$MODE" == "IBM" ]]; then
        [[ "$dir" != *"/ibm/"* ]] && include_dir=false
    fi

    if [[ "$MODE" == "ALL" ]]; then
        include_dir=true
    fi

    # ---- Disable chemistry cases ----
    if [[ "$PELEPHYSICS" == "OFF" ]]; then
        if [[ " ${dirs_reacc[@]} " =~ " ${dir} " ]]; then
            include_dir=false
        fi
    fi

    # ---- Store selected dirs ----
    if $include_dir; then
        FILTERED_DIRS+=("$dir")
    fi

done

echo " selected tests --"
for dir in "${FILTERED_DIRS[@]}"; do
	echo $dir 
done	
echo ""


# Default behavior
STOP_IF_FAIL=FALSE

# Parse arguments
for arg in "$@"; do
    case $arg in
        STOP_IF_FAIL=TRUE)
            STOP_IF_FAIL=TRUE
            ;;
        STOP_IF_FAIL=FALSE)
            STOP_IF_FAIL=FALSE
            ;;
    esac
done


# debug commands
#echo $PELEPHYSICS
#echo $STOP_IF_FAIL
#exit 1

	
# make instruction
MAKECOMP="make -j8"

echo "Checking  Cerisse Examples (be patient) ..."
echo "  "
echo "  "

#header
printf "%-40s\t%-15s\t%-15s\n" "Directory" "Clean" "Build"
printf "%-40s\t%-15s\t%-15s\n" "---------" "-----" "-----"

for dir in "${FILTERED_DIRS[@]}"; do

	#echo $dir

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
			if [[ "$STOP_IF_FAIL" == "TRUE" ]]; then
				echo "===== TPL build failed in $dir ====="
				tail -n 5 build.log
				echo "===================================="
				exit 1
			fi
		fi	    	
	fi 

	#echo -n "Running make in $dir..."
	if $MAKECOMP -C "$dir" > build.log 2>&1; then
		#echo -e "                      ${GREEN}[PASS]${NC}"
		build_result="${GREEN}[PASS]${NC}" 
	else 
		build_result="${RED}[FAIL]${NC}"		
		if [[ "$STOP_IF_FAIL" == "TRUE" ]]; then
			echo "===== compile failed in $dir ====="
			tail -n 5 build.log
			echo "===================================="
			exit 1
		fi
	fi
	# remove compiled files, executables and log files
	rm -rf "${dir:?}/tmp_build_dir"
	rm -rf "${dir:?}/*.ex"
	rm  *.log

	printf "%-40s\t%-15b\t%-15b\n" "$dir" "$clean_result" "$build_result"
done		
printf "%-40s\t%-15s\t%-15s\n" "---------" "-----" "-----"
