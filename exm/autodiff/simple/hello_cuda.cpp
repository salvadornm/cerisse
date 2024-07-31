#include <stdio.h>
#include <stdlib.h>

#ifdef __CUDACC__
# define __CUDA__ __host__ __device__
#else
# define __CUDA__
#endif

#ifdef __CUDACC__
__global__ void kernel(int *a, int *b, int *c)
{
    int tid;

    tid = threadIdx.x + blockIdx.x * blockDim.x;
    c[tid] = a[tid] + b[tid];
    printf(" %d ",tid);
}

int main(void)
{
    int *cuda_a, *cuda_b, *cuda_c, a[100], b[100], c[100];

    printf ("Hello from CUDA\n");

    for (int i = 0; i < 100; i++) {
        a[i] = rand() % 1000;
        b[i] = rand() % 100;
    }
    cudaMalloc((void **)&cuda_a, sizeof(int) * 100);
    cudaMalloc((void **)&cuda_b, sizeof(int) * 100);
    cudaMalloc((void **)&cuda_c, sizeof(int) * 100);
    cudaMemcpy(cuda_a, a, sizeof(int) * 100, cudaMemcpyHostToDevice);
    cudaMemcpy(cuda_b, b, sizeof(int) * 100, cudaMemcpyHostToDevice);
    kernel <<< 10, 10 >>> (cuda_a, cuda_b, cuda_c);
    cudaMemcpy(c, cuda_c, sizeof(int) * 100, cudaMemcpyDeviceToHost);
    cudaFree(cuda_a);
    cudaFree(cuda_b);
    cudaFree(cuda_c);
    return (0);
}
#else
void add(int *a, int *b, int *c)
{
    for (int i = 0; i < 100; i++) {
        c[i] = a[i] + b[i];
    }
}

int main(void)
{
    int a[100], b[100], *c;
    
    printf ("Hello from C++\n");

    c = (int *)malloc(sizeof(int) * 100);
    for (int i = 0; i < 100; i++) {
        a[i] = rand() % 1000;
        b[i] = rand() % 100;
    }
    add(a, b, c);
}
#endif