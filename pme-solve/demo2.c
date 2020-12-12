#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pme-solve.h"
#include "array.h"

static double ic(double x) {
    if (0.0 < x && x < 1.0/2.0) {
        return 1.0/2.0;
    }
    else if (-1.0/2.0 < x && x < 0.0) {
        return -1.0/2.0;
    }
    else {
        return 0;
    }
}

static double bc_L(double t) {
    return 0;
}

static double bc_R(double t) {
    return 0;
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
        .error = 0,
        .exact_sol    = NULL,
        .maple_out    = "pme_demo2.mpl",
        .matlab_out   = "pme_demo2.m",
        .geomview_out = "pme_demo2.gv",
    };
    make_matrix(prob.u, m+1, n+1);
    pme_solve(&prob);
    free_matrix(prob.u);

    return EXIT_SUCCESS;
}
