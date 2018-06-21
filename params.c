#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "params.h"

static void print_usage()
{
    fprintf(stderr,
            "Preconditioned Conjgate Gradient Solver Usage:\n"
            "\t-h: print this message (Defaults Values)\n"
						"\t-a: Precision: Single, Double, Both (Both)\n"
            "\t-x: x dimension (100)\n"
            "\t-y: y dimension (100)\n"
            "\t-z: z dimension (100)\n"
            "\t-M: maximum iteration count (300)\n"
            "\t-r: relative residual tolerance (1e-6)\n"
            "\t-p: preconditioner: SSOR, Jacobi, or id (id)\n"
            "\t-w: Helmhotlz coefficient, omega^2 (1.677e-4)\n"
            "\t-l: Helmholtz coeficient, lambda^2 (1.206e3)\n");
}

static void default_params(solve_param_t* params)
{
    params->n_x       = 100;
    params->n_y       = 100;
    params->n_z       = 100;
    params->maxit   = 300;
    params->rtol    = 1e-6;
    params->ptype   = PC_ID;
    params->omega2   = 0.0001677;
    params->lambda2   = 1206.0;
		params->sp_flag=1;
		params->dp_flag=1;
}

int get_params(int argc, char** argv, solve_param_t* params)
{
    extern char* optarg;
    const char* optstring = "ha:x:y:z:M:r:p:w:l:";
    int c;

    #define get_int_arg(c, field) \
        case c: params->field = atoi(optarg); break
    #define get_flt_arg(c, field) \
        case c: params->field = (float) atof(optarg); break

    default_params(params);
    while ((c = getopt(argc, argv, optstring)) != -1) {
        switch (c) {
        case 'h':
            print_usage();
            return -1;
				case 'a':
						if (strcmp(optarg, "Single") == 0) {
                params->sp_flag = 1;
                params->dp_flag = 0;
            } else if (strcmp(optarg, "Double") == 0) {
                params->sp_flag = 0;
                params->dp_flag = 1;
            } else if (strcmp(optarg, "Both") == 0) {
                params->sp_flag = 1;
								params->dp_flag = 1;
            } else {
                fprintf(stderr, "Unknown precision type: %s\n", optarg);
                return -1;
            }
            break;
        case 'p':
            if (strcmp(optarg, "id") == 0) {
                params->ptype = PC_ID;
            } else if (strcmp(optarg, "Jacobi") == 0) {
                params->ptype = PC_Jacobi;
            } else if (strcmp(optarg, "SSOR") == 0) {
                params->ptype = PC_SSOR;
            } else {
                fprintf(stderr, "Unknown preconditioner type: %s\n", optarg);
                return -1;
            }
            break;	
        get_int_arg('x', n_x);
        get_int_arg('y', n_y);
        get_int_arg('z', n_z);
        get_int_arg('M', maxit);
        get_flt_arg('r', rtol);
        get_flt_arg('w', omega2);
        get_flt_arg('l', lambda2);
        default:
            fprintf(stderr, "Unknown option\n");
            return -1;
        }
    }
    return 0;
}
