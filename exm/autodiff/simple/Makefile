cpu:
	clang++-11 -std=c++11 -I /media/monal/main/code/cerisse/lib/install/clad/include -fplugin=/media/monal/main/code/cerisse/lib/install/clad/lib/clad.so test.cpp 

ad_gpu:
  # -x cu specifies 
	nvcc -std=c++20 -x cu -lcudart_static -ldl -lrt -I /media/monal/main/code/cerisse/lib/install/clad/include -L /media/monal/main/code/cerisse/lib/install/clad/lib/clad.so test_cuda.cpp

hello_gpu:
  # -x cu specifies 
	nvcc -std=c++20 -x cu -lcudart_static -ldl -lrt -I /media/monal/main/code/cerisse/lib/install/clad/include -L /media/monal/main/code/cerisse/lib/install/clad/lib/clad.so hello_cuda.cpp -o hello_cuda.out
# clang++-11 -std=c++11 --cuda-path=/usr/local/cuda-12.3 -lcudart_static -ldl -lrt -pthread --cuda-gpu-arch=sm_80 -i /media/monal/main/code/cerisse/lib/install/clad/include -fplugin=/media/monal/main/code/cerisse/lib/install/clad/lib/clad.so test_cuda.cu

#	 -L/usr/local/cuda/lib64 -lcudart_static \
// RUN: -ldl -lrt -pthread -lm -lstdc++ -oForwardMode.out 2>&1 | FileCheck %s