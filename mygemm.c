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
void dgemm2(const double *A, const double *B, double *C, const int n) {
    int i = 0;
    int j = 0;
    int k = 0;
    for (i = 0; i < n; i += 2){
        for (j = 0; j < n; j += 2){
            register double C_0_0 = C[i * n + j];
            register double C_1_0 = i < (n - 1) ? C[(i + 1) * n + j] : 0;
            register double C_0_1 = j < (n - 1) ? C[i * n + (j + 1)] : 0;
            register double C_1_1 = (i < (n - 1)) && (j < (n - 1))? C[(i + 1) * n + (j + 1)] : 0;
            for (k = 0; k < n; k += 2){
                register double A_0_0 = A[i * n + k];
                register double A_1_0 = i < (n - 1) ? A[(i + 1) * n + k] : 0;
                register double A_0_1 = k < (n - 1) ? A[i * n + (k + 1)] : 0;
                register double A_1_1 = (i < (n - 1)) && (k < (n - 1)) ? A[(i + 1) * n + (k + 1)] : 0;
                register double B_0_0 = B[k * n + j];
                register double B_1_0 = k < (n - 1) ? B[(k + 1) * n + j] : 0;
                register double B_0_1 = j < (n - 1) ? B[k * n + (j + 1)] : 0;
                register double B_1_1 = (k < (n - 1)) && (j < (n - 1)) ? B[(k + 1) * n + (j + 1)] : 0;
                C_0_0 += A_0_0 * B_0_0 + A_0_1 * B_1_0;
                C_1_0 += A_1_0 * B_0_0 + A_1_1 * B_1_0;
                C_0_1 += A_0_0 * B_0_1 + A_0_1 * B_1_1;
                C_1_1 += A_1_0 * B_0_1 + A_1_1 * B_1_1;
            }
            C[i * n + j] = C_0_0;
            if (i < (n - 1)) C[(i + 1) * n + j] = C_1_0;
            if (j < (n - 1)) C[i * n + (j + 1)] = C_0_1;
            if (i < (n - 1) && j < (n - 1)) C[(i + 1) * n + (j + 1)] = C_1_1;
        }
    }
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

void bijk(const double *A, const double *B, double *C, const int n, const int b) {
    int i = 0;
    int j = 0;
    int k = 0;
	for(i=0;i<n;i+=b){
		for(j=0;j<n;j+=b){
			for(k=0;k<n;k+=b){
                int iB = i;
				for(iB=i;iB<i+b && iB<n;iB++){
                    int jB = j;
					for(jB=j;jB<j+b && jB<n;jB++){
						register double res=C[iB*n+jB];
                        int kB = k;
						for(kB=k;kB<k+b && kB<n;kB++){
							res+=A[iB*n+kB]*B[kB*n+jB];
                        }
					    C[iB*n+jB]=res;
					}
                }
			}
        }
    }
}

void jik(const double *A, const double *B, double *C, const int n) {
    int j = 0;
    int i = 0;
    int k = 0;
    for(j=0;j<n;j++){
        for(i=0;i<n;i++){
			register double res=0;
			for(k=0;k<n;k++){
				res+=A[i*n+k]*B[k*n+j];
			}
			C[i*n+j]=res;
		}
    }
}

void bjik(const double *A, const double *B, double *C, const int n, const int b) {
    int j = 0;
    int i = 0;
    int k = 0;
    for(j=0;j<n;j+=b){
		for(i=0;i<n;i+=b){	
			for(k=0;k<n;k+=b){
                int jB = j;
				for(jB=j;jB<j+b && jB<n;jB++){
                    int iB = i;
					for(iB=i;iB<i+b && iB<n;iB++){
						register double res=C[iB*n+jB];
                        int kB = k;
						for(kB=k;kB<k+b && kB<n;kB++){
							res+=A[iB*n+kB]*B[kB*n+jB];
                        }
						C[iB*n+jB]=res;
					}
                }
			}
        }
    }
}

void kij(const double *A, const double *B, double *C, const int n) {
    int j = 0;
    int i = 0;
    int k = 0;
    for(k=0;k<n;k++){
		for(i=0;i<n;i++){
			register double res=A[i*n+k];
			for(j=0;j<n;j++){
				C[i*n+j]+=res*B[k*n+j];
			}
		}
    }
}

void bkij(const double *A, const double *B, double *C, const int n, const int b) {
    int j = 0;
    int i = 0;
    int k = 0;
    for(k=0;k<n;k+=b){
		for(i=0;i<n;i+=b){	
			for(j=0;j<n;j+=b){
                int kB = k;
				for(kB=k;kB<k+b && kB<n;kB++){
                    int iB = i;
					for(iB=i;iB<i+b && iB<n;iB++){
						register double res=A[iB*n+kB];
                        int jB = j;
						for(jB=j;jB<j+b && jB<n;jB++){
							C[iB*n+jB]+=res*B[kB*n+jB];
                        }
					}
                }
			}
        }
    }
}


void ikj(const double *A, const double *B, double *C, const int n) {
    int j = 0;
    int i = 0;
    int k = 0;
    for(k=0;k<n;k++){
		for(i=0;i<n;i++){
			register double res=A[i*n+k];
			for(j=0;j<n;j++){
				C[i*n+j]+=res*B[k*n+j];
			}
		}
    }
}

void bikj(const double *A, const double *B, double *C, const int n, const int b) {
    int j = 0;
    int i = 0;
    int k = 0;
    for(i=0;i<n;i+=b){
		for(k=0;k<n;k+=b){	
			for(j=0;j<n;j+=b){
                int iB = i;
				for(iB=i;iB<i+b && iB<n;iB++){
                    int kB = k;
					for(kB=k;kB<k+b && kB<n;kB++){
						register double res=A[iB*n+kB];
                        int jB = j;
						for(jB=j;jB<j+b && jB<n;jB++){
							C[iB*n+jB]+=res*B[kB*n+jB];
                        }
					}
                }
			}
        }
    }
}

void jki(const double *A, const double *B, double *C, const int n) {
    int j = 0;
    int i = 0;
    int k = 0;
    for(j=0;j<n;j++){
		for(k=0;k<n;k++){
			register double res=B[k*n+j];
			for(i=0;i<n;i++){
				C[i*n+j]+=res*A[i*n+k];
			}
		}
    }
}

void bjki(const double *A, const double *B, double *C, const int n, const int b) {
    int j = 0;
    int i = 0;
    int k = 0;
    for(j=0;j<n;j+=b){
		for(k=0;k<n;k+=b){
			for(i=0;i<n;i+=b){
                int jB = j;
				for(jB=j;jB<j+b && jB<n;jB++){
                    int kB = k;
					for(kB=k;kB<k+b && kB<n;kB++){
						register double res=B[kB*n+jB];
                        int iB = i;
						for(iB=i;iB<i+b && iB<n;iB++){
							C[iB*n+jB]+=res*A[iB*n+kB];
                        }
					}
                }
			}
        }
    }
}

void kji(const double *A, const double *B, double *C, const int n) {
    int j = 0;
    int i = 0;
    int k = 0;
    for(j=0;j<n;j++){
		for(k=0;k<n;k++){
			register double res=B[k*n+j];
			for(i=0;i<n;i++){
				C[i*n+j]+=res*A[i*n+k];
			}
		}
    }
}

void bkji(const double *A, const double *B, double *C, const int n, const int b) {
    int j = 0;
    int i = 0;
    int k = 0;
    for(k=0;k<n;k+=b){
		for(j=0;j<n;j+=b){	
			for(i=0;i<n;i+=b){
                int kB = k;
				for(kB=k;kB<k+b && kB<n;kB++){
                    int jB = j;
					for(jB=j;jB<j+b && jB<n;jB++){
						register double res=B[kB*n+jB];
                        int iB = i;
						for(iB=i;iB<i+b && iB<n;iB++){
							C[iB*n+jB]+=res*A[iB*n+kB];
                        }
					}
                }
			}
        }
    }
}
//Cache Reuse part 3 End 

//Cache Reuse part 4
void optimal(const double* A, const double* B, double *C, const int n, const int b)
{

}