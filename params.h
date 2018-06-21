#ifndef PARAMS_H
#define PARAMS_H

/*@T
 *  * \section{Solver parameters}
 *   * 
 *    * The [[solve_param_t]] structure holds the parameters that
 *     * describe the simulation.  These parameters are filled in
 *      * by the [[get_params]] function.  Details of the parameters
 *       * are described elsewhere in the code.
 *        *@c*/
enum {                /* Types of preconditioners available: */
    PC_ID = 1,        /* 1. Identity                         */
    PC_Jacobi = 2,    /* 2. Block Jacobi                     */
    PC_SSOR = 3    		/* 3. SSOR   			 			              */
};

typedef struct solve_param_t {
    int    n_x;     /* Mesh size in x direction */
    int    n_y;     /* Mesh size in y direction*/
    int    n_z;     /* Mesh size in z direction*/
    int    maxit;   /* Maximum PCG iterations */
    double rtol;    /* Relative residual convergence tolerance */
    int    ptype;   /* Preconditioner type */
    double omega2;  /* Helmholtz Coefficient, omega */
    double lambda2; /* Helmhotlz Coeeficient, lambda */
    int sp_flag;    /* Single Precision  */
		int dp_flag;    /* Double Precision  */
} solve_param_t;

int get_params(int argc, char** argv, solve_param_t* params);

/*@q*/
#endif /* PARAMS_H */
