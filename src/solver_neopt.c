/*
 * Tema 2 ASC
 * 2026 Spring
 */
#include "utils.h"

double* my_solver(int N, double *A, double *B, double *x) {
	
	double *C = calloc(N * N, sizeof(double));
	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){
			for (int k = 0; k < N; k++){
				// fac asa deoarece mi se dau ca vectori 
				C[i*N+j]+= A[k*N+i] * B[k*N+j];
			}
		}
	}
	double *D = calloc(N * N, sizeof(double));
	for (int i = 0; i < N; i++){
		for (int j = i; j < N; j++){
			for (int k = 0; k < N; k++){
				// i <= j 
				D[i*N+j] += C[i*N+k] * C[j*N+k];
			}
		}
	}
	for (int i = 0; i < N; i++){
		for (int j = 0; j < i; j++){
			D[i*N+j] = D[j*N+i];
		}
	}
	double *y = calloc(N, sizeof(double));
	int i = 0;
	while (i < N){
		for (int j = 0; j < N; j++){
			for (int k = j; k < N; k++){
				y[j]+= D[j*N+k] * C[i*N+k];
				if (j != k)
					y[k]+= D[k*N+j] * C[i*N+j];
			}
		}
		i++;
	}
	for (int i = 0; i < N; i++){
		y[i] += x[i];
	}
	free(C);
	free(D);
	return y;
}
