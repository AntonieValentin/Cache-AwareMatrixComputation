/*
 * Tema 2 ASC
 * 2026 Spring
 */
#include "utils.h"

double* my_solver(int N, double *A, double *B, double *x) {
    double *C = calloc(N * N, sizeof(double));
    double *y = calloc(N, sizeof(double));

    // transpus A ca sa nu mai parcurg pe coloana dupa
    double *At = calloc(N * N, sizeof(double));
    for (int i = 0; i < N; i++){
        register double *pointlinieAt = At + i * N;
        for (int j = 0; j < N; j++){
            *(pointlinieAt+j) = A[j*N+i];
        }
    }

    // fac si Bt sa eficientizez
    double *Bt = calloc(N * N, sizeof(double));
    for (int i = 0; i < N; i++){
        register double *pointlinieBt = Bt + i * N;
        for (int j = 0; j < N; j++){
            *(pointlinieBt+j) = B[j*N+i];
        }
    }

    // pun i, j, k
    for (int i = 0; i < N; i++){
        register double *pointlinieC = C + i * N;
        register double *pointlinieA = At + i * N;
        for (int j = 0; j < N; j++){
            register double *pointlinieB = Bt + j * N;
            register double sum = 0;
            for (int k = 0; k < N; k += 40){
                // fac asa deoarece mi se dau ca vectori
                sum += *(pointlinieA+k) * *(pointlinieB+k);
                sum += *(pointlinieA+k+1) * *(pointlinieB+k+1);
                sum += *(pointlinieA+k+2) * *(pointlinieB+k+2);
                sum += *(pointlinieA+k+3) * *(pointlinieB+k+3);
                sum += *(pointlinieA+k+4) * *(pointlinieB+k+4);
                sum += *(pointlinieA+k+5) * *(pointlinieB+k+5);
                sum += *(pointlinieA+k+6) * *(pointlinieB+k+6);
                sum += *(pointlinieA+k+7) * *(pointlinieB+k+7);
                sum += *(pointlinieA+k+8) * *(pointlinieB+k+8);
                sum += *(pointlinieA+k+9) * *(pointlinieB+k+9);
                sum += *(pointlinieA+k+10) * *(pointlinieB+k+10);
                sum += *(pointlinieA+k+11) * *(pointlinieB+k+11);
                sum += *(pointlinieA+k+12) * *(pointlinieB+k+12);
                sum += *(pointlinieA+k+13) * *(pointlinieB+k+13);
                sum += *(pointlinieA+k+14) * *(pointlinieB+k+14);
                sum += *(pointlinieA+k+15) * *(pointlinieB+k+15);
                sum += *(pointlinieA+k+16) * *(pointlinieB+k+16);
                sum += *(pointlinieA+k+17) * *(pointlinieB+k+17);
                sum += *(pointlinieA+k+18) * *(pointlinieB+k+18);
                sum += *(pointlinieA+k+19) * *(pointlinieB+k+19);
                sum += *(pointlinieA+k+20) * *(pointlinieB+k+20);
                sum += *(pointlinieA+k+21) * *(pointlinieB+k+21);
                sum += *(pointlinieA+k+22) * *(pointlinieB+k+22);
                sum += *(pointlinieA+k+23) * *(pointlinieB+k+23);
                sum += *(pointlinieA+k+24) * *(pointlinieB+k+24);
                sum += *(pointlinieA+k+25) * *(pointlinieB+k+25);
                sum += *(pointlinieA+k+26) * *(pointlinieB+k+26);
                sum += *(pointlinieA+k+27) * *(pointlinieB+k+27);
                sum += *(pointlinieA+k+28) * *(pointlinieB+k+28);
                sum += *(pointlinieA+k+29) * *(pointlinieB+k+29);
                sum += *(pointlinieA+k+30) * *(pointlinieB+k+30);
                sum += *(pointlinieA+k+31) * *(pointlinieB+k+31);
                sum += *(pointlinieA+k+32) * *(pointlinieB+k+32);
                sum += *(pointlinieA+k+33) * *(pointlinieB+k+33);
                sum += *(pointlinieA+k+34) * *(pointlinieB+k+34);
                sum += *(pointlinieA+k+35) * *(pointlinieB+k+35);
                sum += *(pointlinieA+k+36) * *(pointlinieB+k+36);
                sum += *(pointlinieA+k+37) * *(pointlinieB+k+37);
                sum += *(pointlinieA+k+38) * *(pointlinieB+k+38);
                sum += *(pointlinieA+k+39) * *(pointlinieB+k+39);
            }
            *(pointlinieC+j) = sum;
        }
    }
    // tin transpusa pentru eficienta
    double *Ct = calloc(N * N, sizeof(double));
    for (int i = 0; i < N; i++){
        register double *pointlinieCt = Ct + i * N;
        for (int j = 0; j < N; j++){
            *(pointlinieCt+j) = C[j*N+i];
        }
    }
    // insumez liniile
    double *sumv = calloc(N, sizeof(double));
    for (int i = 0; i < N; i++){
        register double *pointlinieC = C + i * N;
        for (int j = 0; j < N; j++){
            sumv[j] += *(pointlinieC+j);
        }
    }
    double *temp = calloc(N, sizeof(double));
    for (int i = 0; i < N; i++){
        register double *pointlinieCt = Ct + i * N;
        for (int j = 0; j < N; j++){
            temp[i] += *(pointlinieCt+j) * sumv[j];
        }
    }
    for (int i = 0; i < N; i++){
        register double *pointlinieC = C + i * N;
        for (int j = 0; j < N; j++){
            y[i] += *(pointlinieC+j) * temp[j];
        }
    }
    free(sumv);
    free(temp);
    for (int i = 0; i < N; i++){
        y[i] += x[i];
    }
    free(At);
    free(Bt);
    free(Ct);
    free(C);
    return y;
}