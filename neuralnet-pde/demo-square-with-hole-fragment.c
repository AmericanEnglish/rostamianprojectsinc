#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "neural-nets-pde.h"
#include "array.h"
#include "nelder-mead.h"

static double exact_sol(double x, double y) {
    return (2+x)*(1-pow(x,2))*(1-pow(y,2))*(pow(x,2)+pow(y,2)-1/4);
}
/* The cut-off function for demo-square-with-hole.c is
 *     phi(x,y) = (1-x^2)*(1-y^2)*(x^2 + y^2 - 1/4).
 * The derivatives of phi where calculated symbolically
 * in Maple.  The result was translated to C code also in Maple.
*/
static void my_phi(struct Neural_Net_PDE *nn, double x, double y)
{
	double t1 = x * x;
	double t2 = -t1 + 0.1e1;
	double t4 = y * y;
	double t5 = -t4 + 0.1e1;
	double t6 = t1 + t4 - 0.1e1 / 0.4e1;
	double t7 = t5 * t6;
	double t9 = x * t5;
	double t11 = t2 * t5;
	double t15 = t2 * y;
	double t23 = 0.2e1 * t11;

	nn->phi[0][0] = t2 * t7;
	nn->phi[1][0] = 0.2e1 * (t11 * x - t9 * t6);
	nn->phi[0][1] = 0.2e1 * (t11 * y - t15 * t6);
	nn->phi[2][0] = -0.8e1 * t1 * t5 + t23 - 0.2e1 * t7;
	nn->phi[1][1] = 0.4e1 * (x * y * t6 - t15 * x - t9 * y);
	nn->phi[0][2] = -0.8e1 * t2 * t4 - 0.2e1 * t2 * t6 + t23;
}

/* For the exact solution we pick u(x,t) = (2-x)*phi(x,y),
 * and then we "reverse engineer" the PDE by picking
 * f = u_xx + u_yy.  The calculation of f is not difficult
 * by hand, but I did it in Maple to avoid risk of errors.
*/
static double my_pde(double x, double y, double u,
		double u_x, double u_y,
		double u_xx, double u_xy, double u_yy)
{
	double t1 = x * x;
	double t2 = t1 * t1;
	double t6 = y * y;
	double t16 = t6 * t6;
	double f = 0.2e1 * t2 * x + 0.4e1 * t2
		+ (0.64e2 * t6 - 0.49e2) * t1 * x / 0.2e1
		+ (0.96e2 * t6 - 0.66e2) * t1 / 0.2e1
		+ (0.12e2 * t16 - 0.51e2 * t6 + 0.20e2) * x / 0.2e1
		+ 0.4e1 * t16 - 0.33e2 * t6 + 0.10e2;
	return u_xx + u_yy - f;
}

static void show_usage(char *progname) {
    printf("Usage:\n");
    printf("%s q nu\n", progname);
    printf("q   : number of units in the hidden layer (q >= 1)\n");
    printf("nu  : number of training points (nu >= 1)\n");
}

int main(int argc, char *argv[]) {
    double **training_points = NULL;
    char *endptr;

    if (argc != 3) {
        show_usage(argv[0]);
        return EXIT_FAILURE;
    }

    // number of units in the hidden layer
    int q = strtol(argv[1], &endptr, 10);
    if (*endptr != '\0' || q < 1) {
        show_usage(argv[0]);
        return EXIT_FAILURE;
    }

    // number of training points
    int nu = strtol(argv[2], &endptr, 10);
    if (*endptr != '\0' || nu < 1) {
        show_usage(argv[0]);
        return EXIT_FAILURE;
    }


    struct Neural_Net_PDE nn = {
        .PDE       = my_pde,
        .phi_func  = my_phi,
        .bb_xrange = {-1, 1},
        .bb_yrange = {-1, 1},
        .q         = q,
        .nu        = nu,
        .training_points = training_points,
        .exact_sol       = exact_sol
    };
    printf("q = %d, nu = %d\n", q, nu);
    // nu equally spaced points inside (a,b)
    make_matrix(nn.training_points, nu, 2);
    int count = 0;
    for (int i = 0; i <nu; i++) {
        double r = (double) rand() / RAND_MAX;
        double s = (double) rand() / RAND_MAX;
        double x = (1-r)*nn.bb_xrange[0]+r*nn.bb_xrange[1];
        double y = (1-s)*nn.bb_yrange[0]+s*nn.bb_yrange[1];
        nn.phi_func(&nn, x, y);
        if (nn.phi[0][0] >= 0) {
            nn.training_points[i][0] = x;
            nn.training_points[i][1] = y;
            count++;
        }
    }
    printf("number of training points = %d of %d\n", count, nu);

    Neural_Net_init(&nn);
    struct nelder_mead NM = {
        .f = Neural_Net_residual,
        .n = nn.nweights,
        .s = NULL,
        .x = nn.weights,
        .h = 0.1,
        .tol = 1e-5,
        .maxevals = 1e6,
        .params = &nn,
    };
    printf("weights before training:\n");
    print_vector("%7.3f ", nn.weights, nn.nweights);

    int evalcount = nelder_mead(&NM); // training: minimize the residual
    if (evalcount > NM.maxevals) {
        printf("Nelder-Mead: No convergence after %d function evaluation\n", evalcount);
        return EXIT_FAILURE;
    }
    else {
        printf("Nelder-Mead: Converged after %d function evaluations\n", evalcount);
        printf("Nelder-Mead: Neural network's residual err = %g\n",
                NM.minval);
    }

    printf("weights after training: \n");
    print_vector("%7.3f ", nn.weights, nn.nweights);

    if (nn.exact_sol != NULL) {
        printf("Error versus the PDE's exact solution = %g\n",
                Neural_Net_error_vs_exact(&nn, 50, 50));
                /*Neural_Net_error_vs_exact(&nn, 50, 50));*/
    }
    Neural_Net_plot_with_maple(&nn, 50, 50, "./zz.mpl");
    /*Neural_Net_plot_with_maple(&nn, 50, "/tmp/zz.mpl");*/
    Neural_Net_plot_with_matlab(&nn, 50, 50, "./zz.m");
    /*Neural_Net_plot_with_matlab(&nn, 50, "/tmp/zz.m");*/

    Neural_Net_end(&nn); // end neural network

    free_vector(training_points);
    return EXIT_SUCCESS;
}
