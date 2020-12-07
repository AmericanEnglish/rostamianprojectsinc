#include <stdlib.h>
#include <stdio.h>
#include "heat-solve.h"
#include "array.h"

static void trisolve(int n, const double *a, double *d, const double *c, double *b, double *x) {
    double m;
    for (int i = 1; i < n; i++) {
        m = a[i-1]/d[i-1];
        d[i] -= m*c[i-1];
        b[i] -= m*b[i-1];
    }
    x[n-1] = b[n-1]/d[n-1];
    for (int i = n-2; i >= 0; i--) {
        x[i] = (b[i] - c[i]*x[i+1])/d[i];
    }
}

static double error_vs_exact(struct heat_solve *prob) {
    
}

static void write_plotting_script(struct heat_solve *prob) {
    if (prob->geomview_out != NULL) {
        /*plot_with_geomview(prob);*/
    }
    if (prob->maple_out != NULL) {
        /*plot_with_maple(prob);*/
    }
    if (prob->matlab_out != NULL) {
        /*plot_with_matlab(prob);*/
    }
}

static void heat_solve_implicit(struct heat_solve *prob) {
    int m = prob->m;
    int n = prob->n;
    double **u = prob->u;
    double dx = (prob->b - prob->a) / n;
    double dt = prob->T / m;
    double r = dt/(dx*dx);
    double *b, *d, *c;
    printf("--- Implicit finite difference scheme ---\n");
    printf("r = %g\n", r);

    make_vector(d, n-1);
    make_vector(c, n-2);
    make_vector(b, n+1);

    for (int j =0; j < n-2; j++) {
        c[j] = -r;
    }

    for (int j = 0; j <= n; j++) {
        double x = prob->a + j*dx;
        u[0][j] = prob->ic(x);
    }

    for (int i = 1; i <= m; i++) {
        double t = i*dt;
        u[i][0] = prob->bcL(t);
        u[i][n] = prob->bcR(t);

        for (int j = 1; j <= n-1; j++) {
            b[j] = u[i-1][j];
        }

        b[1]   += r*u[i][0];
        b[n-1] += r*u[i][n];

        for (int k = 0; k < n-1; k++) {
            d[k] = 1+2*r;
        }

        trisolve(n-1, c, d, c, b+1, u[i]+1);
    }

    free_vector(b);
    free_vector(c);
    free_vector(d);
}

void show_usage_and_exit(char *progname) {

}

static void heat_solve_explicit(struct heat_solve *prob) {

}


static void heat_solve_crank_nicolson(struct heat_solve *prob) {

}

static void heat_solve_seidman_sweep(struct heat_solve *prob) {

}

void heat_solve(struct heat_solve *prob) {
    switch (prob->method) {
        case FD_explicit:
            heat_solve_explicit(prob);
            break;
        case FD_implicit:
            heat_solve_implicit(prob);
            break;
        case FD_crank_nicolson:
            heat_solve_crank_nicolson(prob);
            break;
        case FD_seidman_sweep:
            heat_solve_seidman_sweep(prob);
            break;
        default:
            fprintf(stderr, "*** Missing 'method' specifications "
                    "in struct heat_solve\n");
    }

    write_plotting_script(prob);
    if (prob->exact_sol != NULL) {
        prob->error = error_vs_exact(prob);
    }
}
