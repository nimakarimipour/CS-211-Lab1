#include "mygemm.h"

/**
 * 
 * Implement all functions here in this file.
 * Do NOT change input parameters and return type.
 * 
 **/

//Register Reuse part 1
void dgemm0(const double* A, const double* B, double* C, const int n){
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

void dgemm1(const double *A, const double *B, double *C, const int n) {
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
    int i=0;
    int j=0;
    int k=0;
    for (i=0;i<n;i+=2){
        for (j=0;j<n;j+=2){
            register double C00=C[i*n+j];
            register double C10=0;
            if(i<(n-1)){
                C10=C[(i+1)*n+j];
            }
            register double C01=0;
            if(j<(n-1)){
                C01 = C[i*n+(j+1)];
            }
            register double C11=0;
            if((i<(n-1))&&(j<(n-1))){
                C11 = C[(i+1)*n+(j+1)];
            }
            for (k=0;k<n;k += 2){
                register double A00=A[i*n+k];
                register double B00=B[k*n+j];
                register double A10=0; 
                register double A01=0;
                register double B11=0; 
                register double B01=0;
                register double B10=0; 
                register double A11=0; 
                if(k<(n-1)){
                    A01=A[i*n+(k+1)];
                }
                if(i<(n-1)){
                    A10=A[(i+1)*n+k];
                }
                if((i<(n-1))&&(k<(n-1))){
                    A11=A[(i+1)*n+(k+1)];
                }
                if(k<(n-1)){
                    B10=B[(k+1)*n+j];
                }
                if(j<(n-1)){
                    B01=B[k*n+(j+1)];
                }

                if((k<(n-1))&&(j<(n-1))){
                    B11=B[(k+1)*n+(j+1)];
                }
                C00+=A00*B00+A01*B10;
                C01+=A00*B01+A01*B11;
                C10+=A10*B00+A11*B10;
                C11+=A10*B01+A11*B11;
            }
            C[i*n+j]=C00;
            if (i<(n-1))C[(i+1)*n+j]=C10;
            if (j<(n-1))C[i*n+(j+1)]=C01;
            if (i<(n-1)&&j<(n-1))C[(i+1)*n+(j+1)]=C11;
        }
    }
}
//Register Reuse part 2 End

//Register Reuse part 3
void dgemm3(const double *A, const double *B, double *C, const int n) {
    int i = 0;
    int j = 0;
    int k = 0;
    for(i=0;i<n;i+=3){
        for(j=0;j<n;j+=4){
            register double C00=C[i*n+j];
            register double C10=0;
            register double C20=0;
            register double C01=0; 
            register double C11=0; 
            register double C21=0; 
            register double C02=0; 
            register double C12=0;
            register double C22=0;
            register double C03=0;
            register double C13=0;  
            register double C23=0; 
            if(i<(n-1)){
                C10=C[(i+1)*n+j];
            }
            if(i<(n-2)){
                C20=C[(i+2)*n+j];
            }
            if(j<(n-1)){
                C01=C[i*n+(j+1)];
            }
            if(i<(n-1)&&j<(n-1)){
                C11=C[(i+1)*n+(j+1)];
            }
            if(i<(n-2)&&j<(n-1)){
                C21=C[(i+2)*n+(j+1)];
            }
            if(j<(n-2)){
                C02=C[i*n+(j+2)];
            }
            if(i<(n-1)&&j<(n-2)){ 
                C12=C[(i+1)*n+(j+2)];
            }
            if(i<(n-2)&&j<(n-2)){
                C22=C[(i+2)*n+(j+2)];
            }
            if(j<(n-3)){
                C03=C[i*n+(j+3)];
            }
            if(i<(n-1)&&j<(n-3)){
                C13=C[(i+1)*n+(j+3)];
            }
            if(i<(n-2)&&j<(n-3)){
                C23=C[(i+2)*n+(j+3)];
            }
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
    int kb = 0;
    int jb = 0;
    int ib = 0;
	for(i=0;i<n;i+=b){
		for(j=0;j<n;j+=b){
			for(k=0;k<n;k+=b){
                ib = i;
				for(ib=i;ib<i+b && ib<n;ib++){
                    jb = j;
					for(jb=j;jb<j+b && jb<n;jb++){
						register double res=C[ib*n+jb];
                        kb = k;
						for(kb=k;kb<k+b && kb<n;kb++){
							res+=A[ib*n+kb]*B[kb*n+jb];
                        }
					    C[ib*n+jb]=res;
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
    register double ans;
    for(j=0;j<n;j++){
        for(i=0;i<n;i++){
		    ans=0;
			for(k=0;k<n;k++){
				ans+=A[i*n+k]*B[k*n+j];
			}
			C[i*n+j]=ans;
		}
    }
}

void bjik(const double *A, const double *B, double *C, const int n, const int b) {
    int j = 0;
    int i = 0;
    int k = 0;
    int jb = 0;
    int kb = 0;
    int ib = 0;
    for(j=0;j<n;j+=b){
		for(i=0;i<n;i+=b){	
			for(k=0;k<n;k+=b){
                jb = j;
				for(jb=j;jb<j+b && jb<n;jb++){
                    ib = i;
					for(ib=i;ib<i+b && ib<n;ib++){
						register double res=C[ib*n+jb];
                        kb = k;
						for(kb=k;kb<k+b && kb<n;kb++){
							res+=A[ib*n+kb]*B[kb*n+jb];
                        }
						C[ib*n+jb]=res;
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
    int jb = 0;
    int kb = 0;
    int ib = 0;
    for(k=0;k<n;k+=b){
		for(i=0;i<n;i+=b){	
			for(j=0;j<n;j+=b){
                kb = k;
				for(kb=k;kb<k+b && kb<n;kb++){
                    ib = i;
					for(ib=i;ib<i+b && ib<n;ib++){
						register double res=A[ib*n+kb];
                        jb = j;
						for(jb=j;jb<j+b && jb<n;jb++){
							C[ib*n+jb]+=res*B[kb*n+jb];
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
    int kb = 0;
    int jb = 0;
    int ib = 0;
    for(i=0;i<n;i+=b){
		for(k=0;k<n;k+=b){	
			for(j=0;j<n;j+=b){
                ib = i;
				for(ib=i;ib<i+b && ib<n;ib++){
                    kb = k;
					for(kb=k;kb<k+b && kb<n;kb++){
						register double res=A[ib*n+kb];
                        jb = j;
						for(jb=j;jb<j+b && jb<n;jb++){
							C[ib*n+jb]+=res*B[kb*n+jb];
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
    int kb = 0;
    int jb = 0;
    int ib = 0;
    for(j=0;j<n;j+=b){
		for(k=0;k<n;k+=b){
			for(i=0;i<n;i+=b){
                jb = j;
				for(jb=j;jb<j+b && jb<n;jb++){
                    kb = k;
					for(kb=k;kb<k+b && kb<n;kb++){
						register double res=B[kb*n+jb];
                        ib = i;
						for(ib=i;ib<i+b && ib<n;ib++){
							C[ib*n+jb]+=res*A[ib*n+kb];
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
    int kb = 0;
    int jb = 0;
    int ib = 0;
    for(k=0;k<n;k+=b){
		for(j=0;j<n;j+=b){	
			for(i=0;i<n;i+=b){
                kb = k;
				for(kb=k;kb<k+b && kb<n;kb++){
                    jb = j;
					for(jb=j;jb<j+b && jb<n;jb++){
						register double res=B[kb*n+jb];
                        ib = i;
						for(ib=i;ib<i+b && ib<n;ib++){
							C[ib*n+jb]+=res*A[ib*n+kb];
                        }
					}
                }
			}
        }
    }
}
//Cache Reuse part 3 End 

//Cache Reuse part 4
void optimal(const double* A, const double* B, double *C, const int n, const int b){
    int i = 0;
    for (i = 0; i < n; i += b)
    {
        int j = 0;
        for (j = 0; j < n; j += b)
        {
            int k = 0;
            for (k = 0; k < n; k += b)
            {
                int i1 = 0;
                for (i1 = i; i1 < (i + b > n? n : (i + b)); i1 += 3)
                {
                    int j1 = 0;
                    for (j1 = j; j1 < (j + b > n? n : (j + b)); j1 += 3)
                    {
                        register double C_0_0 = C[i1 * n + j1];
                        register double C_1_0 = C[(i1 + 1) * n + j1];
                        register double C_2_0 = C[(i1 + 2) * n + j1];

                        register double C_0_1 = C[i1 * n + (j1 + 1)];
                        register double C_1_1 = C[(i1 + 1) * n + (j1 + 1)];
                        register double C_2_1 = C[(i1 + 2) * n + (j1 + 1)];

                        register double C_0_2 = C[i1 * n + (j1 + 2)];
                        register double C_1_2 = C[(i1 + 1) * n + (j1 + 2)];
                        register double C_2_2 = C[(i1 + 2) * n + (j1 + 2)];

                        int k1 = 0;
                        for (k1 = k; k1 < (k + b > n? n : (k + b)); k1++)
                        {
                            register double A_0_M = A[i1 * n + k1];
                            register double A_1_M = A[(i1 + 1) * n + k1];
                            register double A_2_M = A[(i1 + 2) * n + k1];

                            register double B_M =  B[k1 * n + j1];
                            C_0_0 += A_0_M * B_M;
                            C_1_0 += A_1_M * B_M;
                            C_2_0 += A_2_M * B_M;

                            B_M = B[k1 * n + (j1 + 1)];
                            C_0_1 += A_0_M * B_M;
                            C_1_1 += A_1_M * B_M;
                            C_2_1 += A_2_M * B_M;

                            B_M = B[k1 * n + (j1 + 2)];
                            C_0_2 += A_0_M * B_M;
                            C_1_2 += A_1_M * B_M;
                            C_2_2 += A_2_M * B_M;

                        }
                        C[i1 * n + j1] = C_0_0;
                        C[(i1 + 1) * n + j1] = C_1_0;
                        C[(i1 + 2) * n + j1] = C_2_0;

                        C[i1 * n + (j1 + 1)] = C_0_1;
                        C[(i1 + 1) * n + (j1 + 1)] = C_1_1;
                        C[(i1 + 2) * n + (j1 + 1)] = C_2_1;

                        C[i1 * n + (j1 + 2)] = C_0_2;
                        C[(i1 + 1) * n + (j1 + 2)] = C_1_2;
                        C[(i1 + 2) * n + (j1 + 2)] = C_2_2;
                    
                    }
                }
            }
        }
    }
}