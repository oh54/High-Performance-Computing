#include "ass1_lib.h"
extern "C" {
#include "cblas.h"
}

// Nat
void matmult_nat(int m, int n, int k, double ** A, double ** B, double ** C){
	for(int i = 0; i < m;i++){
		for(int j = 0; j < n;j++){
			C[i][j] = 0;
			for(int l = 0; l < k;l++){
				C[i][j] += A[i][l]*B[l][j];
			}
		}
	}
}

// library implementation through cblas
/*
void matmult_lib(int m, int n, int k, double ** A, double ** B, double ** C){
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1 , A, k, B, n);
}
*/

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

