#include "clad/Differentiator/Differentiator.h"
#include <iostream>

__global__ void add(double *a, double *b, double *c, int n) {
  int idx = threadIdx.x;
  if (idx < n)
    c[idx] = a[idx] + b[idx];
}

double fn1(double i, double j) {
  double a[500] = {};
  double b[500] = {};
  double c[500] = {};
  int n = 500;

  for (int idx=0; idx<500; ++idx) {
    a[idx] = 7;
    b[idx] = 9;
  }

  double *device_a = nullptr;
  double *device_b = nullptr;
  double *device_c = nullptr;

  cudaMalloc(&device_a, n * sizeof(double));
  cudaMalloc(&device_b, n * sizeof(double));
  cudaMalloc(&device_c, n * sizeof(double));

  cudaMemcpy(device_a, a, n * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(device_b, b, n * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(device_c, c, n * sizeof(double), cudaMemcpyHostToDevice);

  add<<<1, 700>>>(device_a, device_b, device_c, n);

  cudaDeviceSynchronize();

  cudaMemcpy(c, device_c, n * sizeof(double), cudaMemcpyDeviceToHost);

  double sum = 0;
  for (int idx=0; idx<n; ++idx)
    sum += c[idx];
  
  return sum * i + 2 * sum * j;
}

int main() {
  // Call clad to generate the derivative of f wrt x.
  auto f_dx = clad::differentiate(fn1, "i");
  // Execute the generated derivative function.
  // std::cout << f_dx.execute(/*x=*/3, /*y=*/4) << std::endl;
  // // Dump the generated derivative code to standard output.
  // f_dx.dump();

  exit(0);
  return 0;
  }