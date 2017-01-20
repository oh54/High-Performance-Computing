#ifndef _ass1_lib_
#define _ass1_lib_

extern "C" {
void matmult_nat(int m, int n, int k, double * A, double * B, double * C);
void matmult_lib(int m, int n, int k, double * A, double * B, double * C);
__host__
void matmult_gpu1(int m, int n, int k, double * A, double * B, double * C);
__host__
void matmult_gpu2(int m, int n, int k, double * A, double * B, double * C);
//__host__
//void matmult_gpu3(int m, int n, int k, int p, double * A, double * B, double * C);


__global__
void cudaSeq(int m, int n, int k, double * A, double * B, double * C);
__global__
void cudaPar(int m, int n, int k, double * A, double * B, double * C);
__global__
void cudaPar2(int m, int n, int k, int p, double * A, double * B, double * C);
__global__
void cudaPar4(int m, int n, int k, int p, double * A, double * B, double * C);
__global__
void cudaSMEM(int m, int n, int k, double * A, double * B, double * C);


void matmult_blk(int m, int n, int k, double ** A, double ** B, double ** C, int bs);
void matmult_nmk(int m, int n, int k, double ** A, double ** B, double ** C);
void matmult_nkm(int m, int n, int k, double ** A, double ** B, double ** C);
void matmult_mnk(int m, int n, int k, double ** A, double ** B, double ** C);
void matmult_mkn(int m, int n, int k, double ** A, double ** B, double ** C);
void matmult_kmn(int m, int n, int k, double ** A, double ** B, double ** C);
void matmult_knm(int m, int n, int k, double ** A, double ** B, double ** C);

__host__ __device__
void printMat(double *A, int m, int n);
}
#endif
