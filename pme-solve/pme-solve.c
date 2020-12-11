#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pme-solve.h"
#include "plot3d.h"
#include "array.h"
#define phi(x) pow((x), 3)

static double error_vs_exact(struct pme_solve *prob) {

}

static double croot(double k) {

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
        v[0] = prob->bcL(t-dt/2);
        // Forward Sweep
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
            // ... supple the missing code
            u[i][j] = 0;
        }
        u[i][0] = prob->bcL(t);
    }
    free_vector(v);

    // Plot
    if (prob->exact_sol != NULL) {
        prob->error = error_vs_exact(prob);
    }
}

