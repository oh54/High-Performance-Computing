#include "stdio.h"

void jacobi_seq(double ** u, double ** uo, double ** f, int N, double delta2);

void gauss_seidel(double ** u, double ** f, int N, double delta2);

void jacobi_omp(double ** u, double ** uo, double ** f, int N, double delta2);
