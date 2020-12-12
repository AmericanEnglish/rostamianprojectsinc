#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pme-solve.h"
#include "array.h"

// Barenblatt's solution
static double barenblatt(
        double x, double t, double m, double c, double delta) {
    double alpha = 1/ (m-1);
    double beta  = 1/ (m+1);
    double gamma = (m-1)/(2*m*(m+1));
    // An easy drop in to help clean up code
    double tpd = pow(t + delta, beta);
    // The piece inside of the max call
    double mx = c - gamma*pow(x/tpd, 2);
    // final value
    double u = tpd*pow(fmax(mx, 0), alpha);
    return u;
}

// pme1: Barenblatt's solution with m = 3, and special choices c and 
//      delta
static double exact_sol(double x, double t) {
    double c     = sqrt(3)/15;
    double delta = 1.0/75;

    return barenblatt(x, t, 3, c, delta);
}

static double ic(double x) {
    return exact_sol(x, 0);
}

static double bc_L(double t) {
    return exact_sol(-1, t);
}

static double bc_R(double t) {
    return exact_sol(1, t);
}

int main(int argc, char*argv[]) {
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
    // Get the command line arguments
    struct pme_solve prob = {
        .a   = -1,
        .b   =  1,
        .T   =  T,
        .n   = n,
        .m   = m,
        .ic  = ic,
        .bcL = bc_L,
        .bcR = bc_R,
        .u   = NULL,
        .exact_sol    = exact_sol,
        .maple_out    = "pme_demo1.mpl",
        .matlab_out   = "pme_demo1.m",
        .geomview_out = "pme_demo1.gv",
    };
    make_matrix(prob.u, m+1, n+1);
    pme_solve(&prob);
    free_matrix(prob.u);

    return EXIT_SUCCESS;
}
