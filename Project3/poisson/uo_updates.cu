__global__ void update_uo_multi_kernel0(double * d0_u, double * d0_uo, int N){
	int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
	int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
	
	d0_uo[i*N + j] = d0_u[i*N + j];
}

__global__ void update_uo_multi_kernel1(double * d1_u, double * d1_uo, int N){
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
	d1_uo[i*N + j] = d1_u[i*N + j];
}
