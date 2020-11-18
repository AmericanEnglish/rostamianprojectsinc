#include <math.h>
#include <stdio.h>

#include "neural-nets-pde.h"
#include "nelder-mead.h"
#include "array.h"

#define PI 4.0*atan(1.0)

static double exact_sol(double x, double y) {
    return 16*x*(1-x)*y*(1-y);
}

static double my_pde(double x, double y, 
        double u_x,               double u_y, 
        double u_xx, double u_xy, double u_yy) {
    return 32*(pow(x, 2)+pow(y, 2)-x-y);
}

static void my_phi(struct Neural_Net_PDE *nn,
        double x, double y) {
    // phi(x,y)    = x(1-x)y(1-y)
    nn->phi[0][0]  = x*(1-x)*y*(1-y);
    // phi_x(x,y)  = (2x-1)y(y-1)
    nn->phi[1][0]  = (2*x-1)*y*(y-1);
    // phi_xx(x,y) = 2y(y-1)
    nn->phi[2][0]  = 2*y*(y-1);
    // phi_y(x,y)  = x(x-1)(2y-1)
    nn->phi[0][1]  = x*(x-1)*(2*y-1);
    // phi_yy(x,y) = 2x(x-1)
    nn->phi[0][2]  = 2*x*(x-1);
    // phi_xy(x,y) = (2x-1)(2y-1)
    nn->phi[1][1]  = (2*x-1)*(2*y-1);
}

static void show_usage(char *progname) {
    printf("Usage:\n");
    printf("%s q nu\n", progname);
    printf("q   : number of units in the hidden layer (q >= 1)\n");
    printf("nu  : number of training points (nu >= 1)\n");
}

int main(int argc, char *argv[]) {
    double a = 0;  // left end of the interval
    double b = PI; // right end of the interval
    double **training_points;
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
        .bb_xrange = {0, 1},
        .bb_yrange = {0, 1},
        .q         = q,
        .nu        = nu,
        .training_points = NULL,
        .exact_sol       = exact_sol
    };

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
        .maxevals = 1e5,
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
    print_vector("%7.3f", nn.weights, nn.nweights);

    if (nn.exact_sol != NULL) {
        printf("Error versus the PDE's exact solution = %g\n",
                Neural_Net_error_vs_exact(&nn, 50, 50));
    }
    Neural_Net_plot_with_maple(&nn, 50, 50, "./zz.mpl");
    /*Neural_Net_plot_with_maple(&nn, 50, "/tmp/zz.mpl");*/
    Neural_Net_plot_with_matlab(&nn, 50, 50, "./zz.m");
    /*Neural_Net_plot_with_matlab(&nn, 50, "/tmp/zz.m");*/

    Neural_Net_end(&nn); // end neural network

    free_vector(training_points);
    return EXIT_SUCCESS;
}
