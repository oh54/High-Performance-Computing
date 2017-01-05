#include "ass1_lib.h"
#include <math.h>
#include <algorithm>
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
void matmult_lib(int m, int n, int k, double ** A, double ** B, double ** C){
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0 , A[0], k, B[0], n, 0.0, C[0], n);
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
int m2=m%bs;
int n2=n%bs;
int k2=k%bs;
int m1=m-m2;
int n1=n-n2;
int k1=k-k2;



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

