#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "precondit.h"
//preconditioners

void pc_identity_double(int N, double* Ax, double* x){
	memcpy(Ax, x, N*sizeof(double));
}

void pc_identity_single(int N, float* Ax, float* x){
  memcpy(Ax, x, N*sizeof(float));
}


void block_jacobi_double(int n_x, int n_y, int n_z, double A_0, double A_z, double* z, double* r){
	#define r(i,j,k) (r[((i)*n_y*n_z)+((j)*n_z)+(k)])
	#define z(i,j,k) (z[((i)*n_y*n_z)+((j)*n_z)+(k)])
	//note these are c' and d' on algo
	double* c = malloc(n_z*sizeof(double));
	double* d = malloc(n_z*sizeof(double));
	for(int i=0;i<n_x;i++){
		for(int j=0;j<n_y;j++){
			//Modify Coefficients
			c[0]=A_z/A_0;
			for(int k=1; k<n_z; k++){
				c[k]=A_z/(A_0-c[k-1]*A_z);
			}
			d[0]=r(i,j,0)/A_0;
			for(int k=1; k<n_z; k++){
				d[k]=(r(i,j,k)-d[k-1]*A_z)/(A_0-c[k-1]*A_z);
			}
			//Back Sub 
			z(i,j,n_z-1)=d[n_z-1];
			for(int k=n_z-2;k>=0;k--){
				z(i,j,k)=d[k]-c[k]*z(i,j,k+1);
			}
		}
	}
}

void block_jacobi_single(int n_x, int n_y, int n_z, float A_0, float A_z, float* z, float* r){
  #define r(i,j,k) (r[((i)*n_y*n_z)+((j)*n_z)+(k)])
  #define z(i,j,k) (z[((i)*n_y*n_z)+((j)*n_z)+(k)])
	//note these are c' and d' on algo
	float* c = malloc(n_z*sizeof(float));
  float* d = malloc(n_z*sizeof(float));
  for(int i=0;i<n_x;i++){
    for(int j=0;j<n_y;j++){
	//Modify Coefficients
	c[0]=A_z/A_0;
      for(int k=1; k<n_z; k++){
        c[k]=A_z/(A_0-c[k-1]*A_z);
      }
      d[0]=r(i,j,0)/A_0;
      for(int k=1; k<n_z; k++){
        d[k]=(r(i,j,k)-d[k-1]*A_z)/(A_0-c[k-1]*A_z);
      }
			//Back Sub
			z(i,j,n_z-1)=d[n_z-1];
      for(int k=n_z-2;k>=0;k--){
        z(i,j,k)=d[k]-c[k]*z(i,j,k+1);
      }
    }
  }
}
