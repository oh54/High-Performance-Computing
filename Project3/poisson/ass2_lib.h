
#include "stdio.h"

__global__ void jacobi_seq_kernel(double * d_u, double * d_uo, double * d_f, int N, double delta2);

__global__ void jacobi_single_kernel(double * d_u, double * d_uo, double * d_f, int N, double delta2);

__global__ void jacobi_multi_kernel(double * d0_u, double * d0_uo, double * d0_f, int N, double delta2);

