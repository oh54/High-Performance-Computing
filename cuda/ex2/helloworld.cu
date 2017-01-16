#include "stdio.h"
#include "stdlib.h"
#include "helper_cuda.h"
__global__ void helloWorldKernelFunc(void){
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int max_i = gridDim.x * blockDim.x;

  if (i == 100){
    int *a = (int*) 0x10000; *a = 0;
  }

  printf("Hello World! I'm thread %i out of %i in block %i. My global thread id is %i out of %i\n", threadIdx.x, blockDim.x, blockIdx.x, i, max_i);
}

int main(int argc, char **argv){
  //Allocate memory space on host and device
  //h_data = malloc(...);
  //cudaMalloc(...);
  //Transfer data from host to device
  //cudaMemcpy(...);
  // Kernel lauch
  helloWorldKernelFunc<<<4, 64>>>();
  checkCudaErrors(cudaDeviceSynchronize());
  // Transfer results from device to host
  //cudaMemcpy(...);
  // Free memory
  //free(h_data);
  //cudaFree(...);
} 
