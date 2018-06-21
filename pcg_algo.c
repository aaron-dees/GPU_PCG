#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "matvec_mul.h"
#include "precondit.h"
//pcg

double dot_double(int n, const double* x, const double* y)
{
	double result = 0.0;
	for (int i = 0; i < n; ++i)
		result += x[i]*y[i];
	return result;
}

double dot_single(int n, const float* x, const float* y)
{
  float result = 0.0;
  for (int i = 0; i < n; ++i)
    result += x[i]*y[i];
  return result;
}

void daxpy(int N, double alpha, double* x, double* y){
	for (int i = 0; i < N; ++i)
		 y[i] += alpha*x[i];
}

void daypx(int N, double alpha, double* x, double* y){
  for (int i = 0; i < N; ++i) 
     y[i] = x[i]+alpha*y[i];
}

void saxpy(int N, float alpha, float* x, float* y){
  for (int i = 0; i < N; ++i) 
     y[i] += alpha*x[i];
}

void saypx(int N, float alpha, float* x, float* y){
  for (int i = 0; i < N; ++i)
     y[i] = x[i]+alpha*y[i];
}


double pcg_double(int N, int n_x, int n_y, int n_z, double* restrict x, const double* restrict b,double A_0, double A_x, double A_y, double A_z, double h_x2, double h_y2, double h_z2, int maxit, double rtol, int PCType){

	double* r = malloc(N*sizeof(double));
	double* z = malloc(N*sizeof(double));
	double* q = malloc(N*sizeof(double));
	double* p = malloc(N*sizeof(double));
	double rho0     = 0;
	double rho      = 0;
	double rho_prev = 0;
	double rtol2 = rtol*rtol;
	int is_converged = 0;
	int step;
	/*Form Residual*/
	//Matrix vector multiplication
	matvec_mul_double(N, n_x, n_y, n_z,A_0,A_x,A_y,A_z,r,x);
	for (int i=0; i<N; ++i){
		r[i]=b[i]-r[i];
	}
	for(step=0; step<maxit&&!is_converged;++step){
		//preconditioner
		if(PCType==1)
			pc_identity_double(N,z,r);
		if(PCType==2)
			block_jacobi_double(n_x,n_y,n_z,A_0,A_z,z,r);
		rho_prev=rho;
		//vector dot
		rho=dot_double(N,r,z);
		if(step==0){
			rho0=rho;
			memcpy(p, z, N*sizeof(double));
		} else {
			double beta = rho/rho_prev;
			daypx(N, beta, z, p);	
		}
		//matrix vector multiplication
		matvec_mul_double(N,n_x,n_y,n_z,A_0,A_x,A_y,A_z,q,p);
		double alpha = rho/dot_double(N, p, q);
		daxpy(N, alpha, p, x);
		daxpy(N, -alpha, q, r);
		//note using rtol^2 here so norm doesnt need to be computed at every step
		is_converged=(rho/rho0<rtol2);
	}
	printf("%d steps, residual reduction %g (%s tol %g)\n",step, sqrt(rho/rho0), is_converged ? "<=" : ">", rtol);
	free(p);
	free(q);
	free(z);
	free(r);
	return rho/rho0;
}

float pcg_float(int N, int n_x, int n_y, int n_z, float* restrict x, const float* restrict b,float A_0, float A_x, float A_y, float A_z, float h_x2, float h_y2, float h_z2, int maxit, float rtol, int PCType){

  float* r = malloc(N*sizeof(float));
  float* z = malloc(N*sizeof(float));
  float* q = malloc(N*sizeof(float));
  float* p = malloc(N*sizeof(float));
  float rho0     = 0;
  float rho      = 0;
  float rho_prev = 0;
  float rtol2 = rtol*rtol;
  int is_converged = 0;
  int step;
	/*Form Residual*/
  //Matrix vector multiplication
	matvec_mul_single(N, n_x, n_y, n_z,A_0,A_x,A_y,A_z,r,x);
  for (int i=0; i<N; ++i){
    r[i]=b[i]-r[i];
	}
  for(step=0; step<maxit&&!is_converged;++step){
		//preconditioner
		if(PCType==1)
      pc_identity_single(N,z,r);
    if(PCType==2)
      block_jacobi_single(n_x,n_y,n_z,A_0,A_z,z,r);
    rho_prev=rho;
		//vector dot
		rho=dot_single(N,r,z);
    if(step==0){
      rho0=rho;
      memcpy(p, z, N*sizeof(float));
    } else {
      float beta = rho/rho_prev;
			saypx(N, beta, z, p);
    }
		//matrix vector multiplication
		matvec_mul_single(N,n_x,n_y,n_z,A_0,A_x,A_y,A_z,q,p);
    float alpha = rho/dot_single(N, p, q);
		saxpy(N, alpha, p, x);
    saxpy(N, -alpha, q, r);
    is_converged=(rho/rho0<rtol2);
  }
  printf("%d steps, residual reduction %g (%s tol %g)\n",step, sqrt(rho/rho0), is_converged ? "<=" : ">", rtol);
  free(p);
  free(q);
  free(z);
  free(r);
  return rho/rho0;
}
