#ifndef MATVEC_MUL_H
#define MATVEC_MUL_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "matvec_mul.h"

void matvec_mul_double(int N,int n_x, int n_y, int n_z,double A_0, double A_x, double A_y, double A_z, double* restrict Ax,double* restrict x);

void matvec_mul_single(int N,int n_x, int n_y, int n_z,float A_0, float A_x, float A_y, float A_z,float* restrict Ax,float* restrict x);

#endif /* MATVEC_MUL_H */
