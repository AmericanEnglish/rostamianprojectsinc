#include <stdio.h>
#include <stdlib.h>
#include "neural-net-ode.h"
#include "nelder-mead.h"
#include "array.h"

static void sigmoid(double x, double *sigma) {
    // tbd
}
static void Neural_Net_phi(struct Neural_Net_ODE *nn, double x) {
    // tbd
}

static void Neural_Net_init(struct Neural_Net_ODE *nn) {
    nn->nweights = 3*nn->q;
    make_vector(nn->weights, nn->nweights);
    for (i = 0; i < n->nweights; i++) {
        nn->weights[i] = (double)rand()/RAND_MAX - 0.5;
    }
}

void Neural_Net_end(struct Neural_Net_ODE *nn) {
    // tbd
}

void Neural_Net_eval(struct Neural_net_ODE *nn, double x) {
    int q = nn->q;
    double *u = nn->weights;
    double *v = nn->weights +   q;
    double *w = nn->weights + 2*q;

    for (int j = 0; j <= 2; j++) {
        nn->N[j] = 0.0;
    }

    for (int i = 0; i < q; i++) {
        double z = u[i] + w[i]*x;
        sigmoid(z, nn->sigma);
        for (int j = 0; j <= 2; j++) {
            nn->N[j] += v[i]*pow(w[i], j) * nn->sigma[j];
        }
    }
}

static double residual_at_x(struct Neural_Net_ODE *nn, double x) {
    // tbd
}

double Neural_Net_residual(double *weights, int nweights, void *params) {
    struct Neural_Net_ODE *nn = params;
    int nu = nn->nu; // the number of training points
    double sum = 0.0;

    for (int i = 0; i < nn->nweights; i++) {
        nn->weights[i] = weights[i];
    }
    for (int i =0; i < nu; i++) {
        double x = nn->training_points[i];
        double r = residual_at_x(nn, x);
        sum += r*r;
    }
    return sum;
}

double Neural_Net_error_vs_exact(struct  Neural_Net_ODE *nn, int n) {
    // tbd
}
