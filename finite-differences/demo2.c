#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "array.h"
#include "heat-solve.h"

#define PI 4.0*atan(1.0)

static double exact_sol(double x, double t) {
    // Demo 2 [-1,1]
    if (t > 0) {
        return erf(x/(2*sqrt(t)));
    }
    else if (t == 0 && x < 0) {
        return -1;
    }
    else if (t == 0 && x > 0) {
        return 1;
    }
    else { // if (t == 0 && x == 0) {
        return 0;
    }
}

// Calculate initial and boundary conditions from the exact solution
static double ic (double x ){
    return exact_sol(x, 0);
}

static double bc_L(double t) {
    return exact_sol(-1, t);
}

static double bc_R(double t) {
    return exact_sol(1, t);
}

int main(int argc, char *argv[]) {
    char *endptr;
    if (argc != 4) {
        show_usage_and_exit(argv[0]);
    }

    double T = strtod(argv[1], &endptr);
    if (*endptr != '\0' || T <= 0.0) {
        show_usage_and_exit(argv[0]);
    }
    
    int n = strtol(argv[2], &endptr, 10);
    if (*endptr != '\0' || n < 1) {
        show_usage_and_exit(argv[0]);
    }

    int m = strtol(argv[3], &endptr, 10);
    if (*endptr != '\0' || m < 1) {
        show_usage_and_exit(argv[0]);
    }
    
    struct heat_solve prob = {
        .a   = -1,
        .b   =  1,
        .T   =  T,
        .n   =  n,
        .m   =  m,
        .ic  = ic,

        .bcL          = bc_L,
        .bcR          = bc_R,
        .method       = FD_undefined,
        .exact_sol    = exact_sol,
        .maple_out    = NULL,
        .matlab_out   = NULL,
        .geomview_out = NULL,
    };

    make_matrix(prob.u, m+1, n+1);

    /* Solve the problem through for different algorithms.
     * Note: matlab does not allow a hyphen in filenames
     * so we use underscores for all our output filenames.
     */

    putchar('\n');
    prob.method       = FD_explicit;
    prob.maple_out    = "prob2_explicit.mpl";
    prob.matlab_out   = "prob2_explicit.m";
    prob.geomview_out = "prob2_implicit.gv";
    heat_solve(&prob);
    if (prob.exact_sol != NULL) {
        printf("absolute error = %g\n", prob.error);
    }
    putchar('\n');

    prob.method       = FD_implicit;
    prob.maple_out    = "prob2_implicit.mpl";
    prob.matlab_out   = "prob2_implicit.m";
    prob.geomview_out = "prob2_implicit.gv";
    heat_solve(&prob);
    if (prob.exact_sol != NULL) {
        printf("absolute error = %g\n", prob.error);
    }
    putchar('\n');

    prob.method       = FD_crank_nicolson;
    prob.maple_out    = "prob2_crank_nicolson.mpl";
    prob.matlab_out   = "prob2_crank_nicolson.m";
    prob.geomview_out = "prob2_crank_nicolson.gv";
    heat_solve(&prob);
    if (prob.exact_sol != NULL) {
        printf("absolute error = %g\n", prob.error);
    }
    putchar('\n');

    prob.method       = FD_seidman_sweep;
    prob.maple_out    = "prob2_seidman_sweep.mpl";
    prob.matlab_out   = "prob2_seidman_sweep.m";
    prob.geomview_out = "prob2_seidman_sweep.gv";
    heat_solve(&prob);
    if (prob.exact_sol != NULL) {
        printf("absolute error = %g\n", prob.error);
    }
    putchar('\n');

    free_matrix(prob.u);

    return EXIT_SUCCESS;
}



