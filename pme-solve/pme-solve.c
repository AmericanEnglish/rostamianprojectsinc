#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pme-solve.h"
#include "plot3d.h"
#include "array.h"
#define phi(x) pow((x), 3)

void show_usage_and_exit(char *progname) {
    printf("Usage: %s T n m\n", progname);
    printf("   T : solve over 0 <= t <= T\n");
    printf("   n : number of subintervals in the x direction\n");
    printf("   m : number of subintervals in the t direction\n");
    exit(EXIT_FAILURE);
}

static double error_vs_exact(struct pme_solve *prob) {
    double error = 0;
    double diff  = 0;
    double t = 0;
    double x = 0;
    int m = prob->m;
    int n = prob->n;
    double dt = prob->T / m;
    double dx = (prob->b - prob->a) / n;
    for (int i = 0; i <= n; i++) {
        t = i*dt;
        for (int j = 0; j <= m; j++) {
            double x = prob->a + j*dx;
            diff = fabs(prob->u[i][j] - prob->exact_sol(x, t));
            if (diff > error) {
                error = diff;
            }
        }
    }
    return error;
}

static double croot(double k) {
    double gamma = 108*k+12*sqrt(12+81*pow(k, 2));
    gamma = pow(gamma, 1.0/3);
    return gamma;
}

void pme_solve(struct pme_solve *prob) {
    int n = prob->n;
    int m = prob->m;
    double **u = prob->u;
    double dx = (prob->b - prob->a) / n;
    double dt = prob->T/m;
    double r = dt/(2*dx*dx);
    double s = sqrt(r);
    double *v;

    make_vector(v, n+1);
    printf("r = %g\n", r);

    for (int j = 0; j <= n; j++) {
        double x = prob->a + j*dx;
        u[0][j] = prob->ic(x);
    }

    for (int i = 1; i <= m; i++) {
        double t = i*dt;
        // Forward Sweep
        v[0] = prob->bcL(t-dt/2);
        for (int j = 1; j <= n-1; j++) {
            double c =
                 r*phi(v[j-1])    +       u[i-1][j]
                -r*phi(u[i-1][j]) + r*phi(u[i-1][j+1]);
            v[j] = croot(s*c)/s;
        }
        v[n] = prob->bcR(t-dt/2); // not used

        // Backward Sweep
        u[i][n] = prob->bcR(t);
        for (int j=n-1; j >= 1; j--) {
            double c = 
                r*phi(v[j-1]) +       v[j]
               -r*phi(v[j])   + r*phi(u[i][j+1]);
            u[i][j] = croot(s*c)/s;
        }
        u[i][0] = prob->bcL(t);
    }
    free_vector(v);
    // Plot
    struct plot3d plotter = {
	.a=prob->a,		// a < x < b
	.b=prob->b,		// a < x < b
	.T=prob->T,		// 0 < t < T
	.n=prob->n,			// number of x subintervals
	.m=prob->m,			// number of t subintervals
	.u=prob->u,		// (m+1)x(n+1) data array
	.maple_out    = prob->maple_out,	// output file for maple 
	.matlab_out   = prob->matlab_out,	// output file for matlab 
	.geomview_out = prob->geomview_out	// output file for geomview 
    };
    plot3d(&plotter);
    // Plot
    if (prob->exact_sol != NULL) {
        prob->error = error_vs_exact(prob);
    }
}

