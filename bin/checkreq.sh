
echo " >> CERISSE checking basic requirements: "

failures=0

REQ_PROGRAMS=("python3" "cmake" "make")

# Loop through each program

echo "--checking Installation ... "

for PROG in "${REQ_PROGRAMS[@]}"; do

	if ! command -v $PROG &> /dev/null; then
    		echo " [ERROR] $PROG is NOT installed."
			((failures++))
	else
    		echo " [OK] $PROG is installed "
	fi
done


# 
echo "--checking Python version ... "
pythonvers=9
python3 -c "import sys; sys.version_info < (3,$pythonvers) and sys.exit(1)" \
    && echo " [OK] >  3.$pythonvers " || echo " [WARNING] current version < recommended 3.$pythonvers"

echo "--checking cmake version ... "
CMAKE_VERSION=$(cmake --version | head -n 1 | awk '{print $3}' | cut -d. -f1,2)
cmakevers=3.2

if (( $(echo "$CMAKE_VERSION > $cmakevers" | bc -l) )); then
    echo " [OK] >  $cmakevers "
else
    echo " [WARNING] current version ($CMAKE_VERSION) < recommended $cmakevers "
fi

echo "--checking GNU makefile version ... "
MAKE_VERSION=$(make --version | head -n 1 | awk '{print $3}' | cut -d. -f1,2)
makevers=3.8

if (( $(echo "$MAKE_VERSION > $makevers" | bc -l) )); then
    echo " [OK] >  $makevers "
else
    echo " [WARNING] current version ($MAKE_VERSION) < recommended $makevers "
fi
##### 


echo " Looking  for C++ compilers ..."

REQ_COMPILERS=("c++" "g++" "icpc" "g++-14" "g++-11") 
INS_COMPILERS=(0,0,0,0,0)
VER_COMPILERS=(0,0,0,0,0)
TYPE_COMPILER=(0,0,0,0,0)
SANE_COMPILER=(0,0,0,0,0)
OPENMP_COMPILER=(0,0,0,0,0)
SOURCE_FILE="auxcpp/com20.cpp"
STANDARD="c++20"
OpenMPFLAG="-fopenmp"
SOURCEOP_FILE="auxcpp/com20openmp.cpp"

i_compiler=0
installed_compiler=0

for COMP in "${REQ_COMPILERS[@]}"; do
	# check compiler installed
    if  command -v $COMP &> /dev/null; then
		COMP_VERSION=$($COMP --version | head -n 1 | awk '{print $4}' | cut -d. -f1,2)
		INS_COMPILERS[installed_compiler]=$COMP
		VER_COMPILERS[installed_compiler]=$COMP_VERSION
        # check if Apple/ GCC
		compiler_version=$($COMP --version)
		if echo "$compiler_version" | grep -q "Apple"; then
			TYPE_COMPILER[installed_compiler]="Clang"
		elif echo "$compiler_version" | grep -q "GCC"; then
			TYPE_COMPILER[installed_compiler]="GCC"
		else
			TYPE_COMPILER[installed_compiler]="Unknown"
		fi
       
        #check if compiler is sane	
		$COMP -std=$STANDARD "$SOURCE_FILE"  &> /dev/null
		compile_status=$?
		if [ $compile_status -eq 0 ]; then
			SANE_COMPILER[installed_compiler]="OK"
			rm a.out
		else
			SANE_COMPILER[installed_compiler]="FAIL"
		fi

        #check if OpenMP installed and sane	
		$COMP -std=$STANDARD $OpenMPFLAG "$SOURCEOP_FILE"  &> /dev/null
		compile_status=$?
		if [ $compile_status -eq 0 ]; then
			OPENMP_COMPILER[installed_compiler]="OK"
			rm a.out
		else
			OPENMP_COMPILER[installed_compiler]="FAIL"
		fi

        ########
		((installed_compiler++))	   
    fi
	((i_compiler++))
done

echo " Compilers installed: $installed_compiler"
# print
printf "%-10s %-10s %-10s %-10s %-10s\n" "Compiler" "Type" "Version" "Working" "OpenMP"
echo "-------------------------------------------------"
for (( i=0; i<$installed_compiler; i++ )) do
    printf "%-10s %-10s %-10s %-10s %-10s\n" ${INS_COMPILERS[i]}  ${TYPE_COMPILER[i]}  ${VER_COMPILERS[i]} ${SANE_COMPILER[i]} ${OPENMP_COMPILER[i]}
#	echo  "  ${INS_COMPILERS[i]}  Type: ${TYPE_COMPILER[i]} Version: ${VER_COMPILERS[i]}  Working: ${SANE_COMPILER[i]} OpenMP: ${OPENMP_COMPILER[i]} "
done
echo "-------------------------------------------------"


echo " Looking  for MPIC++ compilers ..."

MPI_REQ_COMPILERS=("mpic++" "mpicxx" "mpicc") 
MPI_INS_COMPILERS=(0,0,0,0,0)
MPI_VER_COMPILERS=(0,0,0,0,0)
MPI_TYPE_COMPILER=(0,0,0,0,0)
MPI_SANE_COMPILER=(0,0,0,0,0)
SOURCE_FILE="com20mpi.cpp"

mpi_compiler=0
mpi_installed_compiler=0
mpi_success=0

for COMP in "${MPI_REQ_COMPILERS[@]}"; do
	# check compiler installed
    if  command -v $COMP &> /dev/null; then
		COMP_VERSION=$($COMP --version | head -n 1 | awk '{print $4}' | cut -d. -f1,2)
		MPI_INS_COMPILERS[mpi_installed_compiler]=$COMP
		MPI_VER_COMPILERS[mpi_installed_compiler]=$COMP_VERSION

        # check if Apple/ GCC
		compiler_version=$($COMP --version)

		if echo "$compiler_version" | grep -q "Apple"; then
			MPI_TYPE_COMPILER[mpi_installed_compiler]="Clang"
		elif echo "$compiler_version" | grep -q "GCC"; then
			MPI_TYPE_COMPILER[mpi_installed_compiler]="GCC"
		else
			MPI_TYPE_COMPILER[mpi_installed_compiler]="Unknown"
		fi
       
        #check if compiler is sane	
		$COMP -std=$STANDARD "$SOURCE_FILE"  &> /dev/null
		compile_status=$?
		if [ $compile_status -eq 0 ]; then
			MPI_SANE_COMPILER[mpi_installed_compiler]="OK"
			((mpi_success++))
			rm a.out
		else
			MPI_SANE_COMPILER[mpi_installed_compiler]="FAIL"
		fi

        ########
		((mpi_installed_compiler++))	   
    fi
	((mpi_compiler++))
done

echo " MPIC++ compilers installed: $mpi_installed_compiler"
printf "%-10s %-10s %-10s %-10s\n" "MPICompiler" "Type" "Version" "Working"
echo "------------------------------------------"
for (( i=0; i<$mpi_installed_compiler; i++ )) do
	printf "%-10s %-10s %-10s %-10s\n" ${MPI_INS_COMPILERS[i]} ${MPI_TYPE_COMPILER[i]} ${MPI_VER_COMPILERS[i]} ${MPI_SANE_COMPILER[i]}
done
echo "------------------------------------------"

if [ "$mpi_success" -eq 0 ]; then
	((failures++))
	echo "[ERROR] no working MPI compiler found "  
fi


echo "--checking for CUDA compiler ..."
CUDA_COMP="nvcc"


if  command -v $CUDA_COMP &> /dev/null; then
	echo " CUDA installed "
else
	echo " no CUDA compiler found "
fi


##### 
echo " <<<<<<<<<<<<<<<  "
if [ $failures -gt 0 ]; then
    echo "-- Errors found : Installation will not work (check errors)"
else
   echo " -- Sucess: Installation may work"
fi
