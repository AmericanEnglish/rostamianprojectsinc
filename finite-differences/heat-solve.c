#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "heat-solve.h"
#include "array.h"
#include "plot3d.h"


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
    for (int i = 0; i <= prob->m; i++) {
        double t = i*dt;
        for (int j = 0; j < prob->n; j++) {
            double x = prob->a+j*dx;
            // math something?
            diff = fabs(prob->u[i][j] - prob->exact_sol(x, t));
            if  (diff > err) {
                err = diff;
            }

        }
    }
    return err;
}

/*static void write_plotting_script(struct heat_solve *prob) {*/
/*}*/

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
    exit(EXIT_FAILURE);
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
        // Compute Ax+b
        for (int j = 1; j <= n-1; j++) {
            u[i][j] = r*u[i-1][j-1] 
                    + (1-2*r) * u[i-1][j] 
                    + r*u[i-1][j+1];
        }
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
        b[1]   += r*u[i][0];
        b[n-1] += r*u[i][n];
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

    make_vector(v, n);
    
    // Initial condition
    for (int j = 0; j <= n; j++) {
        double x = prob->a + j*dx;
        u[0][j] = prob->ic(x);
    }
    for (int i = 1; i <= m; i++) {
        double t = i*dt;
        // Compute v via a forward sweep for t-dt
        v[0] = prob->bcL(t-dt/2);
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
            u[i][j] /= (1+rp);
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

    /*write_plotting_script(prob);*/
    struct plot3d p = {
        .a = prob->a,		// a < x < b
        .b = prob->b,		// a < x < b
        .T = prob->T,		// 0 < t < T
        .n = prob->n,			// number of x subintervals
        .m = prob->m,			// number of t subintervals
        .u = prob->u,		// (m+1)x(n+1) data array
        .maple_out = prob->maple_out,	// output file for maple graphics 
        .matlab_out = prob->matlab_out,	// output file for matlab graphics 
        .geomview_out = prob->geomview_out	// output file for geomview graphics 
    };
    plot3d(&p);
    if (prob->exact_sol != NULL) {
        prob->error = error_vs_exact(prob);
    }
}
