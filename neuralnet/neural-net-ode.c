#include <stdio.h>
#include <stdlib.h>
#include "neural-net-ode.h"
#include "nelder-mead.h"
#include "array.h"

static void sigmoid(double x, double *sigma) {
    // sigmoid
    sigma[0] = 1.0/(1+pow(M_E, -x));
    // sigmoid'
    sigma[1] = pow(M_E, -x) / pow((pow(M_E, -x) + 1), 2);
    // sigmoid''
    sigma[2] = 2*pow(M_E, -2*x) / pow(pow(M_E, -x)+1, 3)
        - sigma[1];
}
static void Neural_Net_phi(struct Neural_Net_ODE *nn, double x) {
    // phi -- expanded form
    nn->phi[0] = -x*x+x*(nn.a+nn.b)-nn.a*nn.b;
    // phi'
    nn->phi[1] = -2*x+nn.a+nn.b;
    // phi''
    nn->phi[2] = -2;
}

static void Neural_Net_init(struct Neural_Net_ODE *nn) {
    nn->nweights = 3*nn->q;
    make_vector(nn->weights, nn->nweights);
    for (i = 0; i < n->nweights; i++) {
        nn->weights[i] = (double)rand()/RAND_MAX - 0.5;
    }
}

void Neural_Net_end(struct Neural_Net_ODE *nn) {
    // Free allocated memory
    // free_vector(nn->training_points);
    free_vectors(nn->weights);
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
    Neural_Net_phi();  // phi
    Neural_Net_eval(); // N
    nn->ODE;           // ODE

    double R = 0.0;
    // Get N
    Neural_Net_eval(nn, x);
    // Get Phi
    Neural_Net_phi(nn, x);
    return nn->ODE(x,                   
        nn->phi[0]*nn->N[0], 
        nn->phi[1]*nn->N[0] + nn->phi[0]*nn->N[1],
        nn->phi[2]*nn->N[0] + 2*nn->phi[1]*nn->N[1] + phi[0]*N[2]
        );
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
    // 
    double m = 0;
    double d = 0;
    double x;
    for (int i = 0; i < n; i++) {
        x = nn->training_points[i]
        Neural_Net_eval(nn, x);
        Neural_Net_phi(nn, x);
        // Phi(x_i) * N(x_i) - u(x_i)
        d = n->phi[0] * n->N[0] - nn->exact_sol(x);
        d = abs(d);
        if (m > d) {
            m = d;
        }
    }   
    return m;
}
