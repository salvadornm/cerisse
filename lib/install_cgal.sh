cd ./cgal
mkdir build
cd ./build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/home/monal/code/cfd/amrsolver/lib/install/cgal ..
make install


#https://www.boost.org/doc/libs/1_66_0/more/getting_started/unix-variants.html#easy-build-and-install
#cd ./boost
#./bootstrap.sh --prefix=
#./b2 install
