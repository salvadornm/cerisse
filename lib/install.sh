AMREXVERSION=25.09
#AMREXVERSION=23.11
PELEPVERSION=25.04
#PELEPVERSION=23.03
BOOSTVERSION=1.81.0
SUNDIALS_VER=7.5.0
#SUNDIALS_VER=6.5.0
CGALVERSION=6.0.1
#CGALVERSION=5.6.1
#dir= $PWD
case  $1 in
	git )
		echo " installing git versions .. "
		gh repo clone AMReX-Combustion/PelePhysics
		;;
	safe )
    rm -rf amrex
		echo " installing AMREX release version .." $AMREXVERSION
		wget https://github.com/AMReX-Codes/amrex/archive/refs/tags/$AMREXVERSION.zip
		unzip $AMREXVERSION.zip
    mv amrex-$AMREXVERSION amrex
		rm $AMREXVERSION.zip
    rm -rf PelePhysics
    echo " installing PelePhysics release version .." $PELEPVERSION
    wget https://github.com/AMReX-Combustion/PelePhysics/archive/refs/tags/v$PELEPVERSION.zip
    unzip v$PELEPVERSION.zip
    mv PelePhysics-$PELEPVERSION PelePhysics
    rm v$PELEPVERSION.zip
    ;;
	amrex)
		rm -rf amrex
		echo " installing AMREX release version .." $AMREXVERSION
		wget https://github.com/AMReX-Codes/amrex/archive/refs/tags/$AMREXVERSION.zip
		unzip $AMREXVERSION.zip
    mv amrex-$AMREXVERSION amrex
		rm $AMREXVERSION.zip	
		;;
	pelephys)
		rm -rf PelePhysics
    echo " installing PelePhysics release version .." $PELEPVERSION
    wget https://github.com/AMReX-Combustion/PelePhysics/archive/refs/tags/v$PELEPVERSION.zip        
    unzip v$PELEPVERSION.zip
    mv PelePhysics-$PELEPVERSION PelePhysics 
		rm v$PELEPVERSION.zip
		;;
  sundials)
    rm -rf sundials    
    echo " downloading SUNDIALS release version .." $SUNDIAL_VER
    wget https://github.com/LLNL/sundials/archive/refs/tags/v$SUNDIALS_VER.tar.gz
    tar -xzf v$SUNDIALS_VER.tar.gz
    mv sundials-$SUNDIALS_VER sundials
    rm v$SUNDIALS_VER.tar.gz
    echo " NOTE: install is done case-by-case by doing: make TPL" 
    ;;  
  cgal)
    case $2 in
      download)        
        rm -rf cgal
        echo " downloading CGAL release version .." $CGALVERSION
        wget https://github.com/CGAL/cgal/archive/refs/tags/v$CGALVERSION.zip
        unzip v$CGALVERSION.zip
        mv cgal-$CGALVERSION cgal
        rm v$CGALVERSION.zip
        echo -e "\\033[1;32m CGAL download \\033[0m"
        rm -rf boost
        echo " downloading boost release version .." $BOOSTVERSION
        wget https://github.com/boostorg/boost/releases/download/boost-$BOOSTVERSION/boost-$BOOSTVERSION.zip
        unzip boost-$BOOSTVERSION.zip
        mv boost-$BOOSTVERSION boost
        rm boost-$BOOSTVERSION.zip
        echo -e "\\033[1;32m Boost download \\033[0m"
        ;;
      install)
        mkdir -p install
        cd $PWD/boost
        ./bootstrap.sh --prefix=$PWD/../install/boost --with-toolset=gcc
        ./b2 install     
        cd ../cgal
        mkdir -p build
        cd build
        cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$PWD/../../install/cgal ..
        make install
        echo -e "\\033[1;32m CGAL installed \\033[0m"
		    ;;
      *)
		esac
	  ;;
  autodiff)
    git submodule update --init --recursive autodiff
    rm -rf build/autodiff
    mkdir -p ./build/autodiff
    cd build/autodiff
    cmake ../../autodiff/ -DCMAKE_INSTALL_PREFIX=../../install/autodiff
    ;;
  clad)
    git submodule update --init --recursive clad
    mkdir -p build/clad
    mkdir -p install/clad
    cd build/clad
    cmake ../../clad/ -DClang_DIR=/usr/lib/llvm-11 -DLLVM_DIR=/usr/lib/llvm-11 -DCMAKE_INSTALL_PREFIX=../../install/clad -DLLVM_EXTERNAL_LIT="$(which lit)"
    make && make install
    ;;
  *)
    echo " no option selected [git/safe/amrex/pelephys/sundials/cgal]"
    echo "Options: (with default versions 2025)"
    echo "  git           Install using git clone latest AMREX+PelePhysics"    
    echo "  safe          Install using release versions of AMREX+PelePhysics"
    echo "  amrex         Install AMREX release version: $AMREXVERSION "
    echo "  pelephys      Install PelePhysics release version: $PELEPVERSION "
    echo "  sundials      Download SUNDIALS release version: $SUNDIALS_VER "
    echo "  cgal download Download CGAL and Boost release versions: $CGALVERSION and $BOOSTVERSION "
    echo "  cgal install  Install CGAL and Boost in install/ directory "
    exit
esac
echo -e "\\033[1;32m  Installation done \\033[0m"
