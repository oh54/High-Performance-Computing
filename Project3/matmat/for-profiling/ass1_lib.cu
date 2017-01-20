#include "stdio.h"
#include "ass1_lib.h"
#include <math.h>
#include <algorithm>
#include "helper_cuda.h"
#include "cublas_v2.h"
extern "C" {
#include "cblas.h"
//}

// Nat
void matmult_nat(int m, int n, int k, double * A, double * B, double * C){
	for(int i = 0; i < m;i++){
		for(int j = 0; j < n;j++){
			C[i*n + j] = 0;
			for(int l = 0; l < k;l++){
				C[i*n + j] += A[i*k + l]*B[l*n + j];
			}
		}
	}
}

// library implementation through cblas
void matmult_lib(int m, int n, int k, double * A, double * B, double * C){
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0 , A, k, B, n, 0.0, C, n);
//printMat(C,m,n);
}

// Sequential on the GPU
__host__
void matmult_gpu1(int m, int n, int k, double * A, double * B, double * C){

	double *d_A, *d_B, *d_C;
	cudaMalloc(&d_A,n*m*sizeof(double));
	cudaMalloc(&d_B,k*m*sizeof(double));
	cudaMalloc(&d_C,n*k*sizeof(double));

 	cudaMemcpy(d_A,A,n*m*sizeof(double), cudaMemcpyHostToDevice);
 	cudaMemcpy(d_B,B,k*m*sizeof(double), cudaMemcpyHostToDevice);
	
	cudaSeq<<<1,1>>>(m,n,k,d_A,d_B,d_C);
	cudaMemcpy(C,d_C,n*m*sizeof(double), cudaMemcpyDeviceToHost);

	cudaFree(d_A);
	cudaFree(d_B);
	cudaFree(d_C);

}

// Naive GPU - 1 thread per element in C
__host__
void matmult_gpu2(int m, int n, int k, double * A, double * B, double * C){
	int K = 16;
	int gridx = ceil(n*1.0/K);
	int gridy = ceil(m*1.0/K);
	double *d_A, *d_B, *d_C;
	cudaMalloc(&d_A,k*m*sizeof(double));
	cudaMalloc(&d_B,k*n*sizeof(double));
	cudaMalloc(&d_C,n*m*sizeof(double));

 	cudaMemcpy(d_A,A,k*m*sizeof(double), cudaMemcpyHostToDevice);
 	cudaMemcpy(d_B,B,k*n*sizeof(double), cudaMemcpyHostToDevice);
	cudaPar<<<dim3(gridx,gridy),dim3(K,K)>>>(m,n,k,d_A,d_B,d_C);
	cudaMemcpy(C,d_C,n*m*sizeof(double), cudaMemcpyDeviceToHost);
/*	#ifndef __print
	#define __print 5
	printMat(C,m,n);
	#endif
*/
	cudaFree(d_A);
	cudaFree(d_B);
	cudaFree(d_C);

}

// Each thread does neighbouring fields
__host__
void matmult_gpu3(int m, int n, int k, double * A, double * B, double * C){
	int K = 16;
	int p = 2;
	int gridx = ceil(n*1.0/K);
	int gridy = ceil(m*1.0/K/p);
	double *d_A, *d_B, *d_C;

	checkCudaErrors(cudaMalloc(&d_A,k*m*sizeof(double)));
	checkCudaErrors(cudaMalloc(&d_B,n*k*sizeof(double)));
	checkCudaErrors(cudaMalloc(&d_C,n*m*sizeof(double)));

 	checkCudaErrors(cudaMemcpy(d_A,A,k*m*sizeof(double), cudaMemcpyHostToDevice));
 	checkCudaErrors(cudaMemcpy(d_B,B,k*n*sizeof(double), cudaMemcpyHostToDevice));
	cudaPar2<<<dim3(gridx,gridy),dim3(K,K)>>>(m,n,k,p,d_A,d_B,d_C);
//	checkCudaErrors(cudaDeviceSynchronize());
	cudaMemcpy(C,d_C,n*m*sizeof(double), cudaMemcpyDeviceToHost);
//	printMat(C,m,n);

	cudaFree(d_A);
	cudaFree(d_B);
	cudaFree(d_C);
}


void matmult_gpu4(int m, int n, int k, double * A, double * B, double * C){
	int K = 16;
	int p = 4;
	int gridx = ceil(n*1.0/K);
	int gridy = ceil(m*1.0/K/p);
	double *d_A, *d_B, *d_C;

	checkCudaErrors(cudaMalloc(&d_A,k*m*sizeof(double)));
	checkCudaErrors(cudaMalloc(&d_B,n*k*sizeof(double)));
	checkCudaErrors(cudaMalloc(&d_C,n*m*sizeof(double)));

 	checkCudaErrors(cudaMemcpy(d_A,A,k*m*sizeof(double), cudaMemcpyHostToDevice));
 	checkCudaErrors(cudaMemcpy(d_B,B,k*n*sizeof(double), cudaMemcpyHostToDevice));
	cudaPar4<<<dim3(gridx,gridy),dim3(K,K)>>>(m,n,k,p,d_A,d_B,d_C);
//	checkCudaErrors(cudaDeviceSynchronize());
	cudaMemcpy(C,d_C,n*m*sizeof(double), cudaMemcpyDeviceToHost);
//	printMat(C,m,n);

	cudaFree(d_A);
	cudaFree(d_B);
	cudaFree(d_C);
}


void matmult_gpu5(int m, int n, int k, double * A, double * B, double * C){
	
//	cudaSetDevice(4);
	int K = 16;
	int gridx = floor(n*1.0/K);
	int gridy = floor(m*1.0/K);
	double *d_A, *d_B, *d_C;

	checkCudaErrors(cudaMalloc(&d_A,k*m*sizeof(double)));
	checkCudaErrors(cudaMalloc(&d_B,n*k*sizeof(double)));
	checkCudaErrors(cudaMalloc(&d_C,n*m*sizeof(double)));

 	checkCudaErrors(cudaMemcpy(d_A,A,k*m*sizeof(double), cudaMemcpyHostToDevice));
 	checkCudaErrors(cudaMemcpy(d_B,B,k*n*sizeof(double), cudaMemcpyHostToDevice));
	cudaSMEM<<<dim3(gridx,gridy),dim3(K,K)>>>(m,n,k,d_A,d_B,d_C);
	cudaDeviceSynchronize();
//	checkCudaErrors(cudaDeviceSynchronize());
	cudaMemcpy(C,d_C,n*m*sizeof(double), cudaMemcpyDeviceToHost);
	//printMat(C,m,n);
	// printf("\n");

	cudaFree(d_A);
	cudaFree(d_B);
	cudaFree(d_C);
}

void matmult_gpulib(int m, int n, int k, double * A, double * B, double * C){
	double alpha = 1.0;
	double beta = 0.0;
	const double *alphap, *betap;
	alphap = &alpha;
	betap = &beta;

	double *d_A, *d_B, *d_C;

	checkCudaErrors(cudaMalloc(&d_A,k*m*sizeof(double)));
	checkCudaErrors(cudaMalloc(&d_B,n*k*sizeof(double)));
	checkCudaErrors(cudaMalloc(&d_C,n*m*sizeof(double)));

 	checkCudaErrors(cudaMemcpy(d_A,A,k*m*sizeof(double), cudaMemcpyHostToDevice));
 	checkCudaErrors(cudaMemcpy(d_B,B,k*n*sizeof(double), cudaMemcpyHostToDevice));

	cublasHandle_t handle;
	cublasCreate(&handle);
//	cublasDGEMM(cublasHandle_t handle, CUBLAS_OP_N, CUBLAS_OP_N,                         m, n, k,                           alpha,                           A, int lda,                           B, int ldb,                           beta,                           C, int ldc)

cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, alphap, d_B, n, d_A, k, betap, d_C, n);

cublasDestroy(handle);
cudaMemcpy(C,d_C,m*n*sizeof(double),cudaMemcpyDeviceToHost);
cudaFree(d_A);
cudaFree(d_B);
cudaFree(d_C);
}



__global__
void cudaSeq(int m, int n, int k, double * A, double * B, double * C){

	for(int i = 0; i < n; i++){
		for(int j = 0;j < m;j++){
			C[i*n + j] = 0;
		}
	}

	for(int i = 0; i < m;i++){
		for(int l = 0; l < k;l++){
			for(int j = 0; j < n;j++){
				C[i*n + j] += A[i*k + l]*B[l*n + j];
			}
		}
	}

}

__global__
void cudaPar(int m, int n, int k, double * A, double * B, double * C){
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	int j = blockIdx.y*blockDim.y + threadIdx.y;

	if(i < n && j < m){
		C[i + j*n] = 0;

		for(int l = 0; l < k;l++){
			C[j*n + i] += A[k*j + l]*B[l*n + i];
		}
	}
}

__global__
void cudaPar2(int m, int n, int k, int p, double * A, double * B, double * C){
//	printf("first in cudaPar2\n");
	const int P = 2;
	double C_r[P]={0.0,0.0};

//	int q = m%p;
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	int j = blockIdx.y*blockDim.y + threadIdx.y;
//	printf("before if\n");

	j = j*p; // allows C[i + (0/1) j*n] indexing
	if(i < n && j < m-p + 1){ // only works if n%p = 0

		for(int pp = 0; pp < p; pp++){
			for(int l = 0; l < k; l++){
				C_r[pp] += A[(j+ pp)*k + l]*B[n*l + i];
			}
			C[(j + pp)*n + i] =  C_r[pp];
		}

	}

	if(i < n && j > m-p && j < m){
		for(int pp = 0; pp < p; pp++){
			for(int l = 0; l < k; l++){
				C_r[pp] += A[(j+ pp)*k + l]*B[n*l + i];
			}
			C[(j + pp)*n + i] =  C_r[pp];
		}


	}

}


__global__
void cudaPar4(int m, int n, int k, int p, double * A, double * B, double * C){
	const int P = 4;
	double C_r[P]={0.0,0.0,0.0,0.0};

//	int q = m%p;
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	int j = blockIdx.y*blockDim.y + threadIdx.y;

	j = j*p; // allows C[i + (0/1) j*n] indexing
	if(i < n && j < m-p + 1){
		for(int pp = 0; pp < p; pp++){
			for(int l = 0; l < k; l++){
				C_r[pp] += A[(j+ pp)*k + l]*B[n*l + i];
			}
			C[(j + pp)*n + i] =  C_r[pp];
		}
	}

	if(i < n && j > m-p && j < m){
		for(int pp = 0; pp < p; pp++){
			for(int l = 0; l < k; l++){
				C_r[pp] += A[(j+ pp)*k + l]*B[n*l + i];
			}
			C[(j + pp)*n + i] =  C_r[pp];
		}
	}

}

__global__
void cudaSMEM(int m, int n, int k, double * A, double * B, double * C){
	
	const int K = 16;	
	int kk = k/K;

	int i = blockIdx.x*blockDim.x + threadIdx.x;
	int j = blockIdx.y*blockDim.y + threadIdx.y;
	int iBlock = threadIdx.x;
	int jBlock = threadIdx.y;


	__shared__ double smemA[K][K];
	__shared__ double smemB[K][K];
	__shared__ double smemC[K][K];

	smemC[jBlock][iBlock] = 0;

	for(int q = 0; q < kk; q++){
		smemB[jBlock][iBlock] = B[i + q*K*n + jBlock*n];
		smemA[jBlock][iBlock] = A[iBlock + j*k + q*K];
		__syncthreads();
		for(int z = 0; z < K; z++){
			smemC[jBlock][iBlock] += smemA[jBlock][z]*smemB[z][iBlock];
		}
		__syncthreads();
	}
	C[i + j*n] = smemC[jBlock][iBlock];
}















// blocking
void matmult_blk(int m, int n, int k, double ** A, double ** B, double ** C, int bs) { 	
 	for(int i2 = 0; i2 < m;i2++){
		for(int j2 = 0; j2 < n;j2++){
			C[i2][j2] = 0;		
		}
 	
	} 
int bsi=bs;
int bsj=bs;
int bsl=bs;




for(int i1 = 0; i1 < m;i1+=bsi){
	if(m-i1 < bs) {bsi=m-i1;
	}
	for(int l1 = 0; l1 < k;l1+=bsl){
		if(k-l1 < bs) {bsl=k-l1;
		}
		for(int j1 = 0; j1 < n;j1+=bsj){
			if(n-j1 < bs) {bsj=n-j1;
			}
			for(int i2 = 0; i2 < bsi; i2++){	
				for(int l2 = 0; l2 < bsl;l2++){	
					for(int j2 = 0; j2 < bsj; j2++){	
							C[i1+i2][j1+j2] += A[i1+i2][l1+l2]*B[l1+l2][j1+j2];
					}
				}
			}
		}
	}
}
/*
for(int i2 = m1; i2 < m; i2++){	
	for(int l2 = n1; l2 < n;l2++){	
		for(int j2 = k1; j2 < k; j2++){	
			C[i2][j2] += A[i2][l2]*B[l2][j2];
			}
		}
	}
*/
} 


// Permutations of kmn
void matmult_kmn(int m, int n, int k, double ** A, double ** B, double ** C){
	for(int i = 0; i < m;i++){
		for(int j = 0; j < n;j++){
			C[i][j] = 0;		
		}
	}

for(int l = 0; l < k;l++){
	for(int i = 0; i < m;i++){
		for(int j = 0; j < n;j++){

				C[i][j] += A[i][l]*B[l][j];
			}
		}
	}
}

void matmult_knm(int m, int n, int k, double ** A, double ** B, double ** C){
	for(int i = 0; i < m;i++){
		for(int j = 0; j < n;j++){
			C[i][j] = 0;		
		}
	}

for(int l = 0; l < k;l++){
	for(int j = 0; j < n;j++){	
		for(int i = 0; i < m;i++){
				C[i][j] += A[i][l]*B[l][j];
			}
		}
	}
}


void matmult_mnk(int m, int n, int k, double ** A, double ** B, double ** C){
	for(int i = 0; i < m;i++){
		for(int j = 0; j < n;j++){
			C[i][j] = 0;		
		}
	}


	for(int i = 0; i < m;i++){
		for(int j = 0; j < n;j++){
			for(int l = 0; l < k;l++){
				C[i][j] += A[i][l]*B[l][j];
			}
		}
	}
}

void matmult_mkn(int m, int n, int k, double ** A, double ** B, double ** C){
	for(int i = 0; i < m;i++){
		for(int j = 0; j < n;j++){
			C[i][j] = 0;		
		}
	}


	for(int i = 0; i < m;i++){
		for(int l = 0; l < k;l++){
			for(int j = 0; j < n;j++){
				C[i][j] += A[i][l]*B[l][j];
			}
		}
	}
}

void matmult_nkm(int m, int n, int k, double ** A, double ** B, double ** C){
	for(int i = 0; i < m;i++){
		for(int j = 0; j < n;j++){
			C[i][j] = 0;		
		}
	}
for(int j = 0; j < n;j++){
	for(int l = 0; l < k;l++){
		for(int i = 0; i < m;i++){
				C[i][j] += A[i][l]*B[l][j];
			}
		}
	}
}

void matmult_nmk(int m, int n, int k, double ** A, double ** B, double ** C){
	for(int i = 0; i < m;i++){
		for(int j = 0; j < n;j++){
			C[i][j] = 0;		
		}
	}

	for(int j = 0; j < n;j++){
		for(int i = 0; i < m;i++){
			for(int l = 0; l < k;l++){
				C[i][j] += A[i][l]*B[l][j];
			}
		}
	}
}

__host__ __device__
void printMat(double *A, int m, int n){
	for(int i = 0; i < m;i++){
		for(int j = 0; j < n;j++){
			printf("%7.0f ",A[i*n + j]);
		}
		printf("\n");
	}	

}

// extern C
}
