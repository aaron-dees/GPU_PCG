//main
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
//#include <math.h>
//#include <assert.h>
#include "params.h"
#include "pcg_algo.h"



void setup_rhs1_double(int n_x, int n_y, int n_z, double* b)
{
    int N = n_x*n_y*n_z;
    memset(b, 0, N*sizeof(double));
    for (int i = 0; i < n_x; ++i) {
        double x = 1.0*(i+1)/(n_x+1);
        for (int j = 0; j < n_y; ++j) {
            double y = 1.0*(i+1)/(n_y+1);
            for (int k = 0; k < n_z; ++k) {
                double z = 1.0*(i+1)/(n_z+1);
                b[(k*n_x+j)*n_y+i] = x*(1-x) * y*(1-y) * z*(1-z);
            }
        }
    }
}

void setup_rhs1_single(int n_x, int n_y, int n_z, float* b)
{
    int N = n_x*n_y*n_z;
    memset(b, 0, N*sizeof(float));
    for (int i = 0; i < n_x; ++i) {
        float x = 1.0*(i+1)/(n_x+1);
        for (int j = 0; j < n_y; ++j) {
            float y = 1.0*(i+1)/(n_y+1);
            for (int k = 0; k < n_z; ++k) {
                float z = 1.0*(i+1)/(n_z+1);
                b[(k*n_x+j)*n_y+i] = x*(1-x) * y*(1-y) * z*(1-z);
            }
        }
    }
}

int main(int argc, char** argv){
	//declare structure to hold parameters and call function
	solve_param_t params;
	if (get_params(argc, argv, &params))
		return -1;
	//set local variables
	int n_x = params.n_x;
  int n_y = params.n_y;
  int n_z = params.n_z;
	int maxit   = params.maxit;
	int PCType = params.ptype;
	int sp_flag=params.sp_flag;
	int dp_flag=params.dp_flag;
	int N = n_x*n_y*n_z;
	if(dp_flag==1){
  	printf("------------------Using Double Precision------------------\n");
		double rtol = params.rtol;
		double omega2 = params.omega2;
		double lambda2 = params.lambda2;
		double h_x2 = 1.0/(n_x*n_x);
		double h_y2 = 1.0/(n_y*n_y);
		double h_z2 = 1.0/(n_z*n_z);
		double A_0 = (2.0*omega2*((1/h_x2)+(1/h_y2)+(lambda2/h_z2))+1.0);
		double A_x = -(omega2)/h_x2;
		double A_y = -(omega2)/h_y2;
		double A_z = -((omega2)*(lambda2))/h_z2;
		printf("Matrix Elements: A_0=%g,A_x=%g,A_y=%g,A_z=%g\n",A_0,A_x,A_y,A_z);
		//allocate space for b, x and r, set eqaul to zero
  	double* b = malloc(N*sizeof(double));
  	double* x = malloc(N*sizeof(double));
  	double* r = malloc(N*sizeof(double));
  	memset(b, 0, N*sizeof(double));
  	memset(x, 0, N*sizeof(double));
  	memset(r, 0, N*sizeof(double));

		//setup RHS
		setup_rhs1_double(n_x,n_y,n_z,b);
		//call pcg code
		struct timeval start,end;
		gettimeofday(&start,NULL);
		pcg_double(N, n_x, n_y, n_z, x, b, A_0,A_x,A_y,A_z,h_x2,h_y2,h_z2, maxit, rtol, PCType);
		gettimeofday(&end,NULL);
		printf("Time taken:%lfs\n", ((end.tv_sec+(double)end.tv_usec/1000000.0)-(start.tv_sec+(double)start.tv_usec/1000000.0)));	

		//Check answer
		matvec_mul_double(N, n_x, n_y, n_z, A_0, A_x, A_y, A_z, r,x);
		double rnorm2=0.0;
		for (int i = 0; i < N; ++i) r[i] = b[i]-r[i];
		for (int i = 0; i < N; ++i) rnorm2 += r[i]*r[i];
		printf("CHECK: resulting rnorm = %g\n", sqrt(rnorm2)); 
		
		free(r);
		free(x);
		free(b);
	}
	if(sp_flag==1){
		printf("------------------Using Single Precision------------------\n");
  	float rtol_s = params.rtol;
		float omega2_s = params.omega2;
  	float lambda2_s = params.lambda2;
  	float h_x2_s = 1.0/(n_x*n_x);
  	float h_y2_s = 1.0/(n_y*n_y);
  	float h_z2_s = 1.0/(n_z*n_z);
  	float A_0_s = (2.0*omega2_s*((1/h_x2_s)+(1/h_y2_s)+(lambda2_s/h_z2_s))+1.0);
  	float A_x_s = -(omega2_s)/h_x2_s;
  	float A_y_s = -(omega2_s)/h_y2_s;
  	float A_z_s = -((omega2_s)*(lambda2_s))/h_z2_s;
		printf("Matrix Elements: A_0=%g,A_x=%g,A_y=%g,A_z=%g\n",A_0_s,A_x_s,A_y_s,A_z_s);
		float* b_s = malloc(N*sizeof(float));
  	float* x_s = malloc(N*sizeof(float));
  	float* r_s = malloc(N*sizeof(float));
  	memset(b_s, 0, N*sizeof(float));
  	memset(x_s, 0, N*sizeof(float));
  	memset(r_s, 0, N*sizeof(float));
		setup_rhs1_single(n_x,n_y,n_z,b_s);
		struct timeval start_s,end_s;
  	gettimeofday(&start_s,NULL);
  	pcg_float(N, n_x, n_y, n_z, x_s, b_s, A_0_s,A_x_s,A_y_s,A_z_s,h_x2_s,h_y2_s,h_z2_s, maxit, rtol_s, PCType);
  	gettimeofday(&end_s,NULL);
  	printf("Time taken:%lfs\n", ((end_s.tv_sec+(double)end_s.tv_usec/1000000.0)-(start_s.tv_sec+(double)start_s.tv_usec/1000000.0)));
		matvec_mul_single(N, n_x, n_y, n_z, A_0_s, A_x_s, A_y_s, A_z_s, r_s,x_s);
  	float rnorm2_s=0.0;
  	for (int i = 0; i < N; ++i) r_s[i] = b_s[i]-r_s[i];
  	for (int i = 0; i < N; ++i) rnorm2_s += r_s[i]*r_s[i];
  	printf("CHECK: resulting rnorm = %g\n", sqrt(rnorm2_s));

  	free(r_s);
  	free(x_s);
  	free(b_s);
  }
	return 0;
}

