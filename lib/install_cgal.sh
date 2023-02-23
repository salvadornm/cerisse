cd ./cgal
mkdir build
cd ./build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/home/monal/code/cfd/amrsolver/lib/install/cgal ..
make install

#cd ./boost
#./bootstrap.sh --prefix=
