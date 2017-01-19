
#include "stdio.h"

__global__ void jacobi_seq_kernel(double * d_u, double * d_uo, double * d_f, int N, double delta2);

__global__ void jacobi_single_kernel(double * d_u, double * d_uo, double * d_f, int N, double delta2);

__global__ void jacobi_multi_kernel0(double * d_u, double * d_uo, double * d_f, int N, double delta2);

__global__ void jacobi_multi_kernel1(double * d_u, double * d_uo, double * d_f, int N, double delta2);
