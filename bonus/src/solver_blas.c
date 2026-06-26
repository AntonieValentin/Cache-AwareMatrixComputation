/*
 * Tema 2 ASC
 * 2024 Spring
 */
#include "utils.h"
#include <cblas.h>

double* my_solver(int N, double *A, double *B, double *x) {
	double *C = calloc(N * N, sizeof(double));
	cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, N, N, N, 1.0, A, N, B, N, 0.0, C, N);
	double *D = calloc(N * N, sizeof(double));
	cblas_dsyrk(CblasRowMajor, CblasUpper, CblasNoTrans, N, N, 1.0, C, N, 0.0, D, N);
	double *y = calloc(N, sizeof(double));
	int i = 0;
	while (i < N){
		cblas_dsymv(CblasRowMajor, CblasUpper, N, 1.0, D, N, &C[i*N], 1, 1.0, y, 1);
		i++;
	}
	cblas_daxpy(N, 1.0, x, 1, y, 1);
	free(C);
	free(D);
	return y;
}
