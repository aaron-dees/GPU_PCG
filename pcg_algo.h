#ifndef PCG_ALGO_H
#define PCG_ALGO_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "matvec_mul.h"
#include "precondit.h"

double pcg_double(int N, int n_x, int n_y, int n_z, double* restrict x, const double* restrict b, double A_0, double A_x, double A_y, double A_z,double h_x2, double h_y2, double h_z2, int maxit, double rtol, int PCType);

double pcg_float(int N, int n_x, int n_y, int n_z, float* restrict x, const float* restrict b, float A_0, float A_x, float A_y, float A_z,float h_x2, float h_y2, float h_z2, int maxit, float rtol, int PCType);

#endif /* PCG_ALGO_H */
