
#include "stdio.h"

__global__ void jacobi_seq_kernel(double * d_u, double * d_uo, double * d_f, int N, double delta2);

__global__ void jacobi_single_kernel(double * d_u, double * d_uo, double * d_f, int N, double delta2);

__global__ void jacobi_multi_kernel0(double * d0_u, double * d0_uo, double * d0_f, int N, double delta2);

__global__ void jacobi_multi_kernel1(double * d1_u, double * d1_uo, double * d1_f, int N, double delta2);

__global__ void update_uo_multi_kernel(double * d0_u, double * d0_uo, int N);

__global__ void update_uo_multi_kernel1(double * d1_u, double * d1_uo, int N);
