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

}

void dgemm1(const double *A, const double *B, double *C, const int n) 
{

}
//Register Reuse part 1 End

//Register Reuse part 2
void dgemm2(const double *A, const double *B, double *C, const int n) 
{

}
//Register Reuse part 2 End

//Register Reuse part 3
void dgemm3(const double *A, const double *B, double *C, const int n) 
{


}
//Register Reuse part 3 End

//Cache Reuse part 3
void ijk(const double *A, const double *B, double *C, const int n){
    int i = 0;
    int j = 0;
    int k = 0;
    for(i=0;i<n;i++){
		for(j=0;j<n;j++)
		{
			register double res=0;
			for(k=0;k<n;k++)
			{
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
    int i = 0;
    int j = 0;
    int k = 0;
    for(j=0;j<n;j++){
		for(i=0;i<n;i++)
		{
			register double res=0;
			for(k=0;k<n;k++)
			{
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