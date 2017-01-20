#include "stdio.h"
#include "stdlib.h"
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include "ass2_lib.h"
#include <omp.h>
#include "helper_cuda.h"


void printMat(double * A, int N){
	int i,j;
	for(i = 0; i < N; i++){
		for(j = 0; j < N; j++){
			printf("%7.0f ",A[i*N + j]);
		}	
		printf("\n");
	}
}

double getMatSum(double * A, int N){
	int i,j;
	double sum = 0.0;
	for(i = 0; i < N; i++){
		for(j = 0; j < N; j++){
			sum += A[i*N + j];
		}	
	}
	return sum;
}



double fnorm_squared(double * u, double * uo, int N){
	int i,j;
	double sum = 0;
	for(i = 1; i <N-1; i++){
		for(j = 1; j<N-1; j++){
			sum += (u[i*N + j]-uo[i*N + j])*(u[i*N + j]-uo[i*N + j]);
		}
	}
	return sum / (N*N);
}

void update_uo(double * u, double * uo, int N){
	int i,j;
	for(i = 0; i<N; i++){
		for(j = 0; j<N; j++){
			uo[i*N + j] = u[i*N + j];
		}
	}
}


void initialize_matrices(double * u, double * uo, double * f, int N, double Nt){

	
	// init loop variables
	int i, j;

	// define uo as zeros\
	uo[x][0] = 0 i.e. outer wall
	for(i = 1; i < N; i++){
		for(j = 1; j < N-1; j++){
			uo[i*N + j] = 0;
		}
	}
	
	// Defining the boundaries to 20
	{
		for(j = 0; j < N; j++) uo[j] = 20;
		for(i = 0; i < N; i++) uo[i*N] = 20;
		for(i = 0; i < N; i++) uo[i*N + N-1] = 20;
	}
		
	//setting u = uo;
	for(i = 0; i<N; i++){
		for(j = 0; j<N; j++){
			u[i*N + j] = uo[i*N + j];
		}
	}
	//printMat(uo,N);

	// defining the f matrix
	for(i = 0; i < N; i++){
		for(j = 0; j < N; j++){
			f[i*N + j] = 0;
		}
	}

	for(i = 4*Nt; i < 5*Nt; i++){
		for(j = 3*Nt; j < 4*Nt; j++){
			f[i*N + j] = 200;
		}
	}
}
__global__
void ChangePointers(double ** p1, double ** p2) 
{
    double * temp = *p1;
    *p1 = *p2;
    *p2 = temp;
}

__global__
void PrintPointers(double * d_u, double * d_uo){
	printf("&D_U: %i\n", d_u);
	printf("&D_UO: %i\n", d_uo);
}



__host__
void doSeq(double * u, double * uo, double * f, int N, double d, int kmax, double delta2, double dd){
	double start = omp_get_wtime(); 
	double *d_u, *d_uo, *d_f;
	int memsize = N*N*sizeof(double);
	cudaMalloc(&d_u, memsize);
	cudaMalloc(&d_uo, memsize);
	cudaMalloc(&d_f, memsize);	
	cudaMemcpy(d_u, u, memsize, cudaMemcpyHostToDevice);
	cudaMemcpy(d_uo, uo, memsize, cudaMemcpyHostToDevice);
	cudaMemcpy(d_f, f, memsize, cudaMemcpyHostToDevice);
	int k = 0;
	double checksum = 1000.0;
	while(k < kmax){
		jacobi_seq_kernel<<<1, 1>>>(d_u, d_uo, d_f, N, delta2);
		double * temp = d_uo;
    		d_uo = d_u;
    		d_u = temp;
		k++;
	}

	cudaMemcpy(uo, d_uo, memsize, cudaMemcpyDeviceToHost);
	printf("%s, ", "CU-SEQ");
	printf("%f, ", omp_get_wtime()-start);
	printf("%i, %.20f, %.0f, %i\n", N, dd, getMatSum(uo, N), k);
	cudaFree(d_u);
	cudaFree(d_uo);
	cudaFree(d_f);
}


__host__
void doSingle(double * u, double * uo, double * f, int N, double d, int kmax, double delta2, double dd){
	double *d_u, *d_uo, *d_f;
	int memsize = N*N*sizeof(double);
	cudaMalloc(&d_u, memsize);
	cudaMalloc(&d_uo, memsize);
	cudaMalloc(&d_f, memsize);	
	cudaMemcpy(d_u, u, memsize, cudaMemcpyHostToDevice);
	cudaMemcpy(d_uo, uo, memsize, cudaMemcpyHostToDevice);
	cudaMemcpy(d_f, f, memsize, cudaMemcpyHostToDevice);
	int k = 0;
	double checksum = 1000.0;	
	cudaSetDevice(6);
	int K = 16;
	int gridx = ceil((N-2)*1.0/(K));
	int gridy = ceil((N-2)*1.0/(K));
	double start = omp_get_wtime(); 
	while(k < kmax){
		jacobi_single_kernel<<<dim3(gridx,gridy),dim3(K,K)>>>(d_u, d_uo, d_f, N, delta2);
		double * temp = d_uo;
    		d_uo = d_u;
    		d_u = temp;
		k++;
	}
	double end = omp_get_wtime(); 
	cudaMemcpy(uo, d_uo, memsize, cudaMemcpyDeviceToHost);

	//printf("MATRIX UO:\n");
	//printMat(uo,N);
	//printf("\n");

	printf("%s, ", "CU-SIN");
	printf("%f, ", end-start);
	printf("%i, %.20f, %i, %.0f\n", N, dd, k, getMatSum(uo, N));
	cudaFree(d_u);
	cudaFree(d_uo);
	cudaFree(d_f);


}


__host__
void doMulti(double * u, double * uo, double * f, int N, double d, int kmax, double delta2, double dd){
	
	double *d0_u, *d0_uo, *d0_f, *d1_u, *d1_uo, *d1_f;
	int memsize = N*N*sizeof(double);
	int Nsize = N*sizeof(double);

	cudaSetDevice(6);
	cudaDeviceEnablePeerAccess(7,0);
	cudaMalloc((void**)&d0_u, memsize/2 + Nsize);
	cudaMalloc((void**)&d0_uo, memsize/2 + Nsize);
	cudaMalloc((void**)&d0_f, memsize/2 + Nsize);
	cudaMemcpy(d0_u, u, memsize/2 + Nsize, cudaMemcpyHostToDevice);
	cudaMemcpy(d0_uo, uo, memsize/2 + Nsize, cudaMemcpyHostToDevice);
	cudaMemcpy(d0_f, f, memsize/2 + Nsize, cudaMemcpyHostToDevice);

	cudaSetDevice(7); 
	cudaDeviceEnablePeerAccess(6,0);
	cudaMalloc((void**)&d1_u, memsize/2 + Nsize);
	cudaMalloc((void**)&d1_uo, memsize/2 + Nsize);
	cudaMalloc((void**)&d1_f, memsize/2 + Nsize);
	cudaMemcpy(d1_u, &u[memsize/2/sizeof(double) -N], memsize/2 + Nsize, cudaMemcpyHostToDevice);
	cudaMemcpy(d1_uo, &uo[memsize/2/sizeof(double) -N], memsize/2 + Nsize, cudaMemcpyHostToDevice);
	cudaMemcpy(d1_f, &f[memsize/2/sizeof(double) -N], memsize/2 + Nsize, cudaMemcpyHostToDevice);
	int k = 0;
	double checksum = 1000.0;	
	int K = 16;
	int gridx = ceil((N-2)*1.0/(K));
 	int gridy = ceil((N-2)*1.0/(K));
	gridy = ceil(gridy*1.0 / 2);
	double start = omp_get_wtime(); 
	while(k < kmax){
		cudaSetDevice(6);
		jacobi_multi_kernel<<<dim3(gridx,gridy),dim3(K,K)>>>(d0_u, d0_uo, d0_f, N, delta2);
		cudaSetDevice(7);
		jacobi_multi_kernel<<<dim3(gridx,gridy),dim3(K,K)>>>(d1_u, d1_uo, d1_f, N, delta2);
		cudaDeviceSynchronize();
		double * temp = d0_uo;
    		d0_uo = d0_u;
    		d0_u = temp;
		double * temp2 = d1_uo;
		d1_uo = d1_u;
    		d1_u = temp2;
		cudaMemcpy(d1_uo, d0_uo+(N-2)/2*N, Nsize, cudaMemcpyDeviceToDevice);
		cudaMemcpy(d0_uo+(N-2)/2*N+N , d1_uo+N, Nsize, cudaMemcpyDeviceToDevice);
		k++;
	}
	double end = omp_get_wtime();
	cudaSetDevice(6);
	cudaMemcpy(uo, d0_uo, memsize/2, cudaMemcpyDeviceToHost);
	cudaSetDevice(7);
	cudaMemcpy(&uo[memsize/2/sizeof(double)], &d1_uo[N], memsize/2, cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();

	//printf("MATRIX UO:\n");
	//printMat(uo,N);
	//printf("\n");

	printf("%s, ", "CU-MUL");
	printf("%f, ", end-start);
	printf("%i, %.20f, %i, %.0f\n", N, dd, k, getMatSum(uo, N));
	cudaFree(d1_u);
	cudaFree(d1_uo);
	cudaFree(d1_f);
	cudaFree(d0_u);
	cudaFree(d0_uo);
	cudaFree(d0_f);

	
}





int main(int argc, char **argv){
	// ./poisson <method type> <NN> <d> <kmax>
	
	int NN;
	double dd;
	int kmax;
	
	sscanf(argv[2] , "%d", &NN);
	sscanf(argv[3] , "%lf", &dd);
	sscanf(argv[4] , "%d", &kmax);

	double d = dd*dd;
	int N = NN + 2;
	double delta = 2.0/N;
	double delta2 = delta*delta; 
	double Nt = N/6.0; 	

	double * u, * uo, * f;

	u = (double*)malloc(N*N*sizeof(double));
	uo = (double*)malloc(N*N*sizeof(double));
	f = (double*)malloc(N*N*sizeof(double));

	initialize_matrices(u,uo,f, N,Nt);
	//printf("MATRIX U:\n");
	//printMat(u,N);
	//printf("\n");
	//printf("MATRIX UO:\n");
	//printMat(uo,N);
	//printf("\n");

	if(strcmp(argv[1], "seq") == 0){
		doSeq(u, uo, f, N, d, kmax, delta2, dd); 
	}

	
	// naive single GPU
	// NN must be multiple of K
	if(strcmp(argv[1], "sin") == 0){
		doSingle(u, uo, f, N, d, kmax, delta2, dd);

	}
	
	
	// naive multi GPU
	// NN must be even multiple of K
	if(strcmp(argv[1], "mul") == 0){
		doMulti(u, uo, f, N, d, kmax, delta2, dd);



	}
	
}















// WEEK2 OPENMP STUFF
/*
	//initialize_matrices(u, uo, f, N, Nt);
	int i, j;
	struct timeval  tv1, tv2;
	double runtime;
	int nruns;

	if(strcmp(argv[1], "jacobi") == 0){

		runtime = 0.0;
		nruns = 0;
		
		while(runtime <= 3.0){
			k = 0;
			checksum = 1000;
			initialize_matrices(u, uo, f, N, Nt);			

			gettimeofday(&tv1, NULL);
			while(checksum > d && k < kmax){
				jacobi_seq(u,uo,f,N,delta2);
				checksum = fnorm_squared(u,uo,N);
				for(i = 0; i<N; i++){
					for(j = 0; j<N; j++){
						uo[i][j] = u[i][j];
					}
				}
				k++;
//				printf("%f \n", checksum);
			}
			gettimeofday(&tv2, NULL);
			runtime += (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 + (double) (tv2.tv_sec - tv1.tv_sec);
			nruns++;
		}
		printf("%s, ", "JAC");
		printf("%f, ", (double) runtime / nruns);
		printf("%i, %.20f, %i, %i\n", N, dd, k, k*nruns);
		
	}

	if(strcmp(argv[1], "gauss") == 0){

		runtime = 0.0;
		nruns = 0;

		while(runtime <= 3.0){
			k = 0;
			checksum = 1000;
			initialize_matrices(u, uo, f, N, Nt);			
			
			gettimeofday(&tv1, NULL);
			while(checksum > d && k < kmax){
				gauss_seidel(u,f,N,delta2);
				checksum = fnorm_squared(u,uo,N);
				for(i = 0; i<N; i++){
					for(j = 0; j<N; j++){
						uo[i][j] = u[i][j];
					}
				}
				k++;
			}
			gettimeofday(&tv2, NULL);
			runtime += (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 + (double) (tv2.tv_sec - tv1.tv_sec);
			nruns++;
		}
		
		printf("%s, ", "G-S");
		printf("%f, ", runtime);
		printf("%i, %.20f, %i, %i\n", N, dd, k, k*nruns);
		
	}
	
	if(strcmp(argv[1], "omp") == 0){
		runtime = 0.0;
	
		k = 0;
		checksum = 1000;
		initialize_matrices(u, uo, f, N, Nt);			
		
		double omp_s = omp_get_wtime();
		while(checksum > d && k < kmax){
			#pragma omp parallel default(none) shared(u,uo,f,N,delta2) private(i,j)
			{
				jacobi_seq(u,uo,f,N,delta2);
			} // end parallel
			checksum = fnorm_squared(u,uo,N);
			for(i = 0; i<N; i++){
				for(j = 0; j<N; j++){
					uo[i][j] = u[i][j];
				}
			} 
			
			k++;
		}
		double omp_time = omp_get_wtime() - omp_s;
		int thread = omp_get_max_threads();
		printf("%s, ", "OMP");
		printf("%f, ", omp_time);
		printf("%i, %.20f, %i, %i\n", N, dd, k, thread);

	}


	if(strcmp(argv[1], "omp2") == 0){
			runtime = 0.0;
	
		k = 0;
		checksum = 1000;
		initialize_matrices(u, uo, f, N, Nt);			
		
		double omp_s = omp_get_wtime();
			#pragma omp parallel default(none) shared	(u,uo,f,N,delta2, checksum, k, d, kmax) private(i,j)
				{
			while(checksum > d && k < kmax){
					
					jacobi_seq(u,uo,f,N,delta2);
					
					// checksum = fnorm_squared(u,uo,N);
					#pragma omp for	private(i,j)  reduction(+:checksum)
					for(i = 1; i <N-1; i++){
						for(j = 1; j<N-1; j++){
							checksum += (u[i][j]-uo[i][j])*(u[i][j]-uo[i][j]);
						}
					}
					
					
					#pragma omp for	private(i,j) 
					for(i = 0; i<N; i++){
						for(j = 0; j<N; j++){
							uo[i][j] = u[i][j];
						}
					} 
					#pragma omp master
					{
					k++;
					checksum=checksum/(N*N);		
					}
					#pragma omp barrier
				} // end while 
				
			} // end parallel
		double omp_time = omp_get_wtime() - omp_s;
		int thread = omp_get_max_threads();
		printf("%s, ", "OMP2");
		printf("%f, ", omp_time);
		printf("%i, %.20f, %i, %i\n", N, dd, k, thread);

	}

	if(strcmp(argv[1], "omp3") == 0){
		runtime = 0.0;
		k = 0;
		checksum = 1000;
		#pragma omp parallel default(none) shared(u, uo, f, N, Nt)
		{
			initialize_matrices(u, uo, f, N, Nt);			
		}
		double omp_s = omp_get_wtime();
			#pragma omp parallel default(none) shared	(u,uo,f,N,delta2, checksum, k, d, kmax) private(i,j)
			{
				while(checksum > d && k < kmax){
					jacobi_seq(u,uo,f,N,delta2);
					// checksum = fnorm_squared(u,uo,N);
					#pragma omp for	private(i,j)  reduction(+:checksum)
					for(i = 1; i <N-1; i++){
						for(j = 1; j<N-1; j++){
							checksum += (u[i][j]-uo[i][j])*(u[i][j]-uo[i][j]);
						}
					}
					#pragma omp for	private(i,j) 
					for(i = 0; i<N; i++){
						for(j = 0; j<N; j++){
							uo[i][j] = u[i][j];
						}
					} 
					#pragma omp master
					{
					k++;
					checksum=checksum/(N*N);
					}
					#pragma omp barrier
				} // end while 
			} // end parallel
	double omp_time = omp_get_wtime() - omp_s;
	int thread = omp_get_max_threads();
	printf("%s, ", "OMP3");
	printf("%f, ", omp_time);
	printf("%i, %.20f, %i, %i\n", N, dd, k, thread);

	}

	

//	printMat(u,N);
	// Save the data

	/*
	The real code should be here. While loop that checks if change from uo to u is small enough to be accepted (solution has converged). Jacobi should be implemented as a sub-routine in a separate function
	*/

	//printf("k is: %i \n",k);

//	dfree_2d(u);
//	dfree_2d(uo);
//	dfree_2d(f);
//}

