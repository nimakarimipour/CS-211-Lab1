#include "mygemm.h"

/**
 * 
 * Implement all functions here in this file.
 * Do NOT change input parameters and return type.
 * 
 **/

//Register Reuse part 1
void dgemm0(const double* A, const double* B, double* C, const int n)
{
    int i = 0;
    int j = 0;
    int k = 0;
    for (i=0; i<n; i++){
        for (j=0; j<n; j++){
            for (k=0; k<n; k++){
                C[i*n+j] += A[i*n+k] * B[k*n+j];
            }
        }
    }
}

void dgemm1(const double *A, const double *B, double *C, const int n) 
{
    int i = 0;
    int j = 0;
    int k = 0;
    for (i=0; i<n; i++){
        for (j=0; j<n; j++) {
            register double r = C[i*n+j]; 
            for (k=0; k<n; k++){
                r += A[i*n+k] * B[k*n+j]; 
            }
            C[i*n+j] = r;
        }
    }
}
//Register Reuse part 1 End

//Register Reuse part 2
void dgemm2(const double *A, const double *B, double *C, const int n) 
{

}
//Register Reuse part 2 End

//Register Reuse part 3
void dgemm3(const double *A, const double *B, double *C, const int n) {
    int i = 0;
    int j = 0;
    int k = 0;
    for (i = 0; i < n; i += 3){
        for (j = 0; j < n; j += 4){
            register double C00 = C[i * n + j];
            register double C10 = i < (n - 1) ? C[(i + 1) * n + j] : 0;
            register double C20 = i < (n - 2) ? C[(i + 2) * n + j] : 0;
            register double C01 = j < (n - 1) ? C[i * n + (j + 1)] : 0;
            register double C11 = i < (n - 1) && j < (n - 1) ? C[(i + 1) * n + (j + 1)] : 0;
            register double C21 = i < (n - 2) && j < (n - 1) ? C[(i + 2) * n + (j + 1)] : 0;
            register double C02 = j < (n - 2) ? C[i * n + (j + 2)] : 0;
            register double C12 = i < (n - 1) && j < (n - 2) ? C[(i + 1) * n + (j + 2)] : 0;
            register double C22 = i < (n - 2) && j < (n - 2) ? C[(i + 2) * n + (j + 2)] : 0;
            register double C03 = j < (n - 3) ? C[i * n + (j + 3)] : 0;
            register double C13 = i < (n - 1) && j < (n - 3) ? C[(i + 1) * n + (j + 3)] : 0;
            register double C23 = i < (n - 2) && j < (n - 3) ? C[(i + 2) * n + (j + 3)] : 0;
            for (k = 0; k < n; k++){
                register double A0M = A[i * n + k];
                register double A1M = i < (n - 1) ? A[(i + 1) * n + k] : 0;
                register double A2M = i < (n - 2) ? A[(i + 2) * n + k] : 0;
                register double BM = B[k * n + j];
                C00 += A0M * BM;
                C10 += A1M * BM;
                C20 += A2M * BM;
                BM = j < (n - 1) ? B[k * n + (j + 1)] : 0;
                C01 += A0M * BM;
                C11 += A1M * BM;
                C21 += A2M * BM;
                BM = j < (n - 2) ? B[k * n + (j + 2)] : 0;
                C02 += A0M * BM;
                C12 += A1M * BM;
                C22 += A2M * BM;
                BM = j < (n - 3) ? B[k * n + (j + 3)] : 0;
                C03 += A0M * BM;
                C13 += A1M * BM;
                C23 += A2M * BM;
            }
            C[i * n + j] = C00;
            if (i < (n - 1)) C[(i + 1) * n + j] = C10;
            if (i < (n - 2)) C[(i + 2) * n + j] = C20;
            if (j < (n - 1)) C[i * n + (j + 1)] = C01;
            if (i < (n - 1) && j < (n - 1)) C[(i + 1) * n + (j + 1)] = C11;
            if (i < (n - 2) && j < (n - 1)) C[(i + 2) * n + (j + 1)] = C21;
            if (j < (n - 2)) C[i * n + (j + 2)] = C02;
            if (i < (n - 1) && j < (n - 2)) C[(i + 1) * n + (j + 2)] = C12;
            if (i < (n - 2) && j < (n - 2)) C[(i + 2) * n + (j + 2)] = C22;
            if (j < (n - 3)) C[i * n + (j + 3)] = C03;
            if (i < (n - 1) && j < (n - 3)) C[(i + 1) * n + (j + 3)] = C13;
            if (i < (n - 2) && j < (n - 3)) C[(i + 2) * n + (j + 3)] = C23;
        }
    }
}
//Register Reuse part 3 End

//Cache Reuse part 3
void ijk(const double *A, const double *B, double *C, const int n){
    int i = 0;
    int j = 0;
    int k = 0;
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            register double res=0;
            for(k=0;k<n;k++){
                res+=A[i*n+k]*B[k*n+j];
            }
            C[i*n+j]=res;
        }       
    }
}

void bijk(const double *A, const double *B, double *C, const int n, const int b) 
{

}

void jik(const double *A, const double *B, double *C, const int n) 
{

}

void bjik(const double *A, const double *B, double *C, const int n, const int b) 
{

}

void kij(const double *A, const double *B, double *C, const int n) 
{

}

void bkij(const double *A, const double *B, double *C, const int n, const int b) 
{

}


void ikj(const double *A, const double *B, double *C, const int n) 
{

}

void bikj(const double *A, const double *B, double *C, const int n, const int b) 
{

}

void jki(const double *A, const double *B, double *C, const int n) 
{

}

void bjki(const double *A, const double *B, double *C, const int n, const int b) 
{

}

void kji(const double *A, const double *B, double *C, const int n) 
{

}

void bkji(const double *A, const double *B, double *C, const int n, const int b) 
{

}
//Cache Reuse part 3 End 

//Cache Reuse part 4
void optimal(const double* A, const double* B, double *C, const int n, const int b)
{

}