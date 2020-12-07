#include <stdlib.h>
#include <stdio.h>
#include "heat-solve.h"
#include "array.h"

#include <math.h>

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
    // iterate over time and space
    double err  = 0;
    double diff;
    double dx = (prob->b - prob->a) / prob->n;
    double dt = (prob->T) / prob->m;
    /* Matrix u is n+1 x m+1 so we should ignore the spatial boundary
     * so we ignore i = 0 and i = n.
     */
    // Ignore the boundary points...
    for (int i = 1; i < prob->m; i++) {
        for (int j = 1; j < prob->n; j++) {
            double x = prob->a+j*dx;
            double t = i*dt;
            // math something?
            diff = fabs(prob->u[i][j] - prob->exact_sol(x, t));
            if  (diff > err) {
                err = diff;
            }

        }
    }
    return err;
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
    // Compute off diagonal
    for (int j =0; j < n-2; j++) {
        c[j] = -r;
    }
    // Initial condition
    for (int j = 0; j <= n; j++) {
        double x = prob->a + j*dx;
        u[0][j] = prob->ic(x);
    }
    // Start time steps
    for (int i = 1; i <= m; i++) {
        // Compute time
        double t = i*dt;
        // Compute boundary condition
        u[i][0] = prob->bcL(t);
        u[i][n] = prob->bcR(t);
        // Fill b in Ax=b
        for (int j = 1; j <= n-1; j++) {
            b[j] = u[i-1][j];
        }
        b[1]   += r*u[i][0];
        b[n-1] += r*u[i][n];
        // Fill the diagonal of A
        for (int k = 0; k < n-1; k++) {
            d[k] = 1+2*r;
        }
        // Solve the system Ax=b
        trisolve(n-1, c, d, c, b+1, u[i]+1);
    }

    free_vector(b);
    free_vector(c);
    free_vector(d);
}

void show_usage_and_exit(char *progname) {
    printf("Usage: %s T n m\n", progname);
    printf("   T : solve over 0 <= T \n");
    printf("   n : grid points a=x[0], ..., b=x[n]\n");
    printf("   m : time slices 0=t[0], ..., T=t[m]\n");

}

static void heat_solve_explicit(struct heat_solve *prob) {
    // Variables Needed
    int m = prob->m;
    int n = prob->n;
    double **u = prob->u;
    double dx = (prob->b - prob->a) / n;
    double dt = prob->T / m;
    double r = dt/(dx*dx);
    printf("--- Explicit finite difference scheme ---\n");
    printf("r = %g\n", r);

    // Initial Condition
    for (int j = 0; j <= n; j++) {
        double x = prob->a + j*dx;
        u[0][j] = prob->ic(x);
    }
    // Time Step
    double t;
    for (int i = 1; i <= m; i++) {
        // Compute time
        t = i*dt;
        // Boundary Conditions
        u[i][0] = prob->bcL(t);
        u[i][n] = prob->bcR(t);
        // Compute Ax
        for (int j = 1; j <= n-1; j++) {
            u[i][j] = r*u[i-1][j-1] 
                    + (1-2*r) * u[i-1][j] 
                    + r*u[i-1][j+1];
        }
        // Compute +b
        u[i][1]   = r*u[i][0];
        u[i][n-1] = r*u[i][n];
    }
}


static void heat_solve_crank_nicolson(struct heat_solve *prob) {
    // Variables Needed
    int m = prob->m;
    int n = prob->n;
    double **u = prob->u;
    double dx = (prob->b - prob->a) / n;
    double dt = prob->T / m;
    double r = dt/(dx*dx);
    double *b, *d, *c;
    printf("--- Crank finite difference scheme ---\n");
    printf("r = %g\n", r);

    make_vector(d, n-1);
    make_vector(c, n-2);
    make_vector(b, n+1);
    // Off Diagonal
    for (int j = 0; j < n-2; j++) {
        c[j] = -r;
    }
    // Initial Condition
    for (int j = 0; j <= n; j++ ) {
        double x = prob->a+j*dx;
        u[0][j] = prob->ic(x);
    }
    // Time Step
    for (int i = 1; i <= m; i++) {
        // Compute time
        double t = i*dt;
        // Boundary Conditions
        u[i][0] = prob->bcL(t);
        u[i][n] = prob->bcR(t);
        // Fill the b in Ax=b
        // In this case we compute the explicit result into b
        for (int j = 1; j <= n-1; j++) {
            // s' = 2(1-r)
            b[j] = r*u[i-1][j-1]
                 + 2*(1-r)*u[i-1][j]
                 + r*u[i-1][j+1];
        }
        // Requires the boundary from the future i guess?
        b[1]   += r*u[i][0]+r*prob->bcL(t+dt);
        b[n-1] += r*u[i][n]+r*prob->bcR(t+dt);
        // Fill diagonal
        for (int k = 0; k < n-1; k++) {
            d[k] = 2*(1+r);
        }
        // Solve the equation for the given timestep
        trisolve(n-1, c, d, c, b+1, u[i]+1);
    }
    free_vector(c);
    free_vector(d);
    free_vector(b);
}

static void heat_solve_seidman_sweep(struct heat_solve *prob) {
    // Variables
    // Variables Needed
    int m = prob->m;
    int n = prob->n;
    double **u = prob->u;
    double dx = (prob->b - prob->a) / n;
    double dt = prob->T / m;
    double r = dt/(dx*dx);
    double rp = r/2;
    double *v;
    printf("--- Seidman finite difference scheme ---\n");
    printf("r = %g\n", r);

    make_vector(v, n-1);
    
    // Initial condition
    for (int j = 0; j <= n; j++) {
        double x = prob->a + j*dx;
        u[0][j] = prob->ic(x);
    }
    for (int i = 1; i <= m; i++) {
        double t = i*dt;
        // Compute v via a forward sweep for t-dt
        for (int j = 1; j <= n-1; j++) {
            // Simplified formula
            v[j] = rp*v[j-1] + (1-rp)*u[i-1][j]+rp*u[i-1][j+1];
            v[j] /= (1+rp);
        }
        // Fill in the boundary points
        u[i][0] = prob->bcL(t);
        u[i][n] = prob->bcR(t);
        // Compute u via the backward sweep for t
        for (int j = n-1; j >= 1; j--) {
            // Simplified formula
            u[i][j] = rp*v[j-1] + (1-rp)*v[j] + rp*u[i][j+1];
        }
    }
    free_vector(v);
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
