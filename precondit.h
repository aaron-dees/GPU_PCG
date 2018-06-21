#ifndef PRECONDIT_H
#define PRECONDIT_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

void pc_identity_double(int N, double* Ax, double* x);
void pc_identity_single(int N, float* Ax, float* x);
void block_jacobi_double(int n_x, int n_y, int n_z, double A_0, double A_z, double* z, double* r);
void block_jacobi_single(int n_x, int n_y, int n_z, float A_0, float A_z, float* z, float* r);
#endif /* PRECONDIT_H */
