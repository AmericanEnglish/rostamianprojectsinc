#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "xmalloc.h"
#include "nelder-mead.h"
#include "neural-nets-ode.h"
#include "array.h"
#define PI 4.0*atan(1.0)

static double exact_sol(double x) {
    return sin(x)+sin(2*x);
}

static double my_ode(double x, double u, double u_x, double u_xx) {
    double t = sin(x) + sin(2*x);
    double f = -sin(x) - 4*sin(2*x) + 1/(1+t*t);
    return u_xx + 1/(1+u*u) - f;
}

static void show_usage(char *progname) {
    printf("Usage:\n");
    printf("%s q nu\n", progname);
    printf("q   : number of units in the hidden layer (q >= 1)\n");
    printf("nu  : number of training points (nu >= 1)\n");
}

int main(int argc, char **argv) {
    double a = 0;  // left end of the interval
    double b = PI; // right end of the interval
    double *training_points;
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

    // nu equally spaced points inside (a,b)
    make_vector(training_points, nu);
    for (int i =0; i < nu; i++) {
        training_points[i] = a + (b - a) / (nu + 1) * (i + 1);
    }
    struct Neural_Net_ODE nn = {
        .a               = a,
        .b               = b,
        .q               = q,
        .nu              = nu,
        .training_points = training_points,
        .ODE             = my_ode,
        .exact_sol       = exact_sol,
    };

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
        printf("Nelder-Mead: No convergence after %d "
                "function evaluation\n", evalcount);
        return EXIT_FAILURE;
    }
    else {
        printf("Nelder-Mead: Converged after %d"
                "function evaluations\n", evalcount);
        printf("Nelder-Mead: Neural network's residual err = %g\n",
                NM.minval);
    }

    printf("weights after training: \n");
    print_vector("%7.3f", nn.weights, nn.nweights);

    if (nn.exact_sol != NULL) {
        printf("Error versus the ODE's exact solution = %g\n",
                Neural_Net_error_vs_exact(&nn, 50));
    }
    Neural_Net_plot_with_maple(&nn, 50, "./zz.mpl");
    /*Neural_Net_plot_with_maple(&nn, 50, "/tmp/zz.mpl");*/
    Neural_Net_plot_with_matlab(&nn, 50, "./zz.m");
    /*Neural_Net_plot_with_matlab(&nn, 50, "/tmp/zz.m");*/

    Neural_Net_end(&nn); // end neural network

    free_vector(training_points);
    return EXIT_SUCCESS;
}
