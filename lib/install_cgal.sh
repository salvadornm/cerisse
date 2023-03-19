cd ../
dir=$PWD

file=$HOME/.bashrc
echo "export AMREX_HOME=$dir/lib/amrex" >> $file

#https://www.boost.org/doc/libs/1_66_0/more/getting_started/unix-variants.html#easy-build-and-install
cd $dir/lib/boost
./bootstrap.sh --prefix=$dir/lib/install/boost
./b2 install

cd $dir/lib/cgal
mkdir build
cd ./build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$dir/lib/install/cgal ..
make install
