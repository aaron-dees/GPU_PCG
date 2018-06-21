#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "matvec_mul.h"
//matrix vector multiplication

void matvec_mul_double(int N,int n_x, int n_y, int n_z,double A_0, double A_x, double A_y, double A_z,double* restrict Ax,double* restrict x){

	#define X(i,j,k) (x[((i)*n_y*n_z)+((j)*n_z)+(k)])
	#define AX(i,j,k) (Ax[((i)*n_y*n_z)+((j)*n_z)+(k)])

	for (int i = 0; i < n_x; ++i) {
		for (int j = 0; j < n_y; ++j) {
			for (int k = 0; k < n_z; ++k) {
				double xx = A_0*X(i,j,k);
				double xn = (i > 0) ? A_x*X(i-1,j,k) : 0;
				double xs = (i < n_x-1) ? A_x*X(i+1,j,k) : 0;
				double xe = (j > 0) ? A_y*X(i,j-1,k) : 0;
				double xw = (j < n_y-1) ? A_y*X(i,j+1,k) : 0;
				double xu = (k > 0) ? A_z*X(i,j,k-1) : 0;
				double xd = (k < n_z-1) ? A_z*X(i,j,k+1) : 0;
				AX(i,j,k) = (xx + xn + xs + xe + xw + xu + xd);
			}
		}
	}
	#undef AX
	#undef X
}

void matvec_mul_single(int N,int n_x, int n_y, int n_z,float A_0, float A_x, float A_y, float A_z,float* restrict Ax,float* restrict x){

  #define X(i,j,k) (x[((i)*n_y*n_z)+((j)*n_z)+(k)])
  #define AX(i,j,k) (Ax[((i)*n_y*n_z)+((j)*n_z)+(k)])

  for (int i = 0; i < n_x; ++i) {
    for (int j = 0; j < n_y; ++j) {
      for (int k = 0; k < n_z; ++k) {
        float xx = A_0*X(i,j,k);
        float xn = (i > 0) ? A_x*X(i-1,j,k) : 0;
        float xs = (i < n_x-1) ? A_x*X(i+1,j,k) : 0;
        float xe = (j > 0) ? A_y*X(i,j-1,k) : 0;
        float xw = (j < n_y-1) ? A_y*X(i,j+1,k) : 0;
        float xu = (k > 0) ? A_z*X(i,j,k-1) : 0;
        float xd = (k < n_z-1) ? A_z*X(i,j,k+1) : 0;
        AX(i,j,k) = (xx + xn + xs + xe + xw + xu + xd);
      }
    }
  }
  #undef AX
  #undef X
}

