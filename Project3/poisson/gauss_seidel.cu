void gauss_seidel(double ** u, double ** f, int N, double delta2){

	int i,j;
	for(i = 1; i < N-1; i++){
		for(j = 1; j < N-1; j++){
			u[i][j] = 0.25*(u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1] + delta2*f[i][j]);
		}
	}
}
