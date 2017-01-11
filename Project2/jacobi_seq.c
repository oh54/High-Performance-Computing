void jacobi_seq(double ** u, double ** uo, double ** f, int N, double delta2){

int i,j;
for(i = 1; i < N-1; i++){
	for(j = 1; j < N-1; j++){
		u[i][j] = 0.25*(uo[i-1][j] + uo[i+1][j] + uo[i][j+1] + uo[i][j-1] + delta2*f[i][j]);
	}
}
}
