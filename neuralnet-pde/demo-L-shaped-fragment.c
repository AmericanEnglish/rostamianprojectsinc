#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "neural-nets-pde.h"
#include "array.h"
#include "nelder-mead.h"

static double exact_sol(double x, double y) {
      double p = (1-pow(x,2))*(1-pow(y,2))
          *(x + y + sqrt(pow(x,2)+pow(y,2)));
      return x*p;
}
/* The cut-off function for demo-L-shaped.c is
 *     phi(x,y) = (1-x^2)*(1-y^2)*(x + y + sqrt(x^2+y^2)).
 * The derivatives of phi where calculated symbolically
 * in Maple.  The result was translated to C code also in Maple.
*/
static void my_phi(struct Neural_Net_PDE *nn, double x, double y)
{
	double t1 = x * x;
	double t2 = t1 - 0.1e1;
	double t4 = y * y;
	double t5 = t4 - 0.1e1;
	double t6 = t1 + t4;
	double t7 = sqrt(t6);
	double t8 = t7 + y + x;
	double t9 = t8 * t5;
	double t11 = t5 * x;
	double t14 = t5 * t2;
	double t15 = 0.1e1 / t7;
	double t17 = t15 * x + 0.1e1;
	double t21 = y * t2;
	double t25 = t15 * y + 0.1e1;
	double t33 = 0.1e1 / t7 / t6;

	nn->phi[0][0] = t9 * t2;
	nn->phi[1][0] = 0.2e1 * t8 * t11 + t17 * t14;
	nn->phi[0][1] = t25 * t14 + 0.2e1 * t8 * t21;
	nn->phi[2][0] = 0.2e1 * t9 + 0.4e1 * t17 * t11 + (-t33 * t1 + t15) * t14;
	nn->phi[1][1] = -y * t33 * x * t14 + 0.4e1 * t8 * y * x + 0.2e1 * t25 * t11 + 0.2e1 * t17 * t21;
	nn->phi[0][2] = 0.2e1 * t8 * t2 + 0.4e1 * t25 * t21 + (-t33 * t4 + t15) * t14;
}

/* For the exact solution we pick u(x,t) = x*phi(x,y),
 * and then we "reverse engineer" the PDE by picking
 * f = u_xx + u_yy.  The calculation of f was done
 * symbolically in Maple.
*/
static double my_pde(double x, double y, double u,
		double u_x, double u_y,
		double u_xx, double u_xy, double u_yy)
{
	double t1 = x * x;
	double t2 = t1 * t1;
	double t6 = y * y;
	double t17 = sqrt(t1 + t6);
	double t22 = t6 * t6;
	double f = 0.2e1 / t17 * (t17 * (t2 + 0.3e1 * y * t1 * x + t1 * (0.6e1 * t6 - 0.7e1) + x * (0.3e1 * t6 * y - 0.6e1 * y) - t6 + 0.1e1) + x * (t2 + t1 * (0.19e2 / 0.2e1 * t6 - 0.15e2 / 0.2e1) + 0.3e1 * t22 - 0.15e2 / 0.2e1 * t6 + 0.3e1 / 0.2e1));
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
