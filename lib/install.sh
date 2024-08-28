AMREXVERSION=23.11
PELEPVERSION=23.03
CGALVERSION=5.6
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
	cgal)
    mkdir -p install
    cd $PWD/boost
    ./bootstrap.sh --prefix=$PWD/../install/boost
    ./b2 install

    cd ../cgal
    mkdir -p build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$PWD/../../install/cgal ..
    make install
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
	echo " no option selected [git/safe/amrex/pelephys/cgal]"
    exit
esac
echo -e "\\033[1;32m  Installation done \\033[0m"
