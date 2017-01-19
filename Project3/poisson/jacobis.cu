#include "helper_cuda.h"
__global__ 
void jacobi_seq_kernel(double * d_u, double * d_uo, double * d_f, int N, double delta2){

	int i,j;
	for(i = 1; i < N-1; i++){
		for(j = 1; j < N-1; j++){
			d_u[i*N + j] = 0.25*(d_uo[(i-1)*N + j] + d_uo[(i+1)*N + j] + d_uo[i*N + j+1] + d_uo[i*N + j-1] + delta2*d_f[i*N + j]);
		}
	}
}

__global__ 
void jacobi_single_kernel(double * d_u, double * d_uo, double * d_f, int N, double delta2){
	int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
	int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
	d_u[i*N + j] = 0.25*(d_uo[(i-1)*N + j] + d_uo[(i+1)*N + j] + d_uo[i*N + j+1] + d_uo[i*N + j-1] + delta2*d_f[i*N + j]);
}

__global__ 
void jacobi_multi_kernel0(double * d_u, double * d_uo, double * d_f, int N, double delta2){
	int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
	int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
	//printf("I AND J 0 (%i, %i)\n", i, j);	
	d_u[i*N + j] = 0.25*(d_uo[(i-1)*N + j] + d_uo[(i+1)*N + j] + d_uo[i*N + j+1] + d_uo[i*N + j-1] + delta2*d_f[i*N + j]);
}

__global__ 
void jacobi_multi_kernel1(double * d_u, double * d_uo, double * d_f, int N, double delta2){
	int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	//printf("I AND J 1 (%i, %i)\n", i, j);	
	d_u[i*N + j] = 0.25*(d_uo[(i-1)*N + j] + d_uo[(i+1)*N + j] + d_uo[i*N + j+1] + d_uo[i*N + j-1] + delta2*d_f[i*N + j]);
}


