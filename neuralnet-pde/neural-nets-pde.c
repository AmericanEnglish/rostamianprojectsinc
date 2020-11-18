#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "neural-nets-pde.h"
#include "array.h"

static void sigmoid(double x, double *sigma) {
    // sigmoid
    sigma[0] = 1.0/(1+exp(-x));
    // sigmoid'
    sigma[1] = exp(-x) / pow((exp(-x) + 1), 2);
    // sigmoid''
    sigma[2] = 2*exp(-2*x) / pow((exp(-x)+ 1), 3)
        - sigma[1];
}

void Neural_Net_init(struct Neural_Net_PDE *nn) {
    nn->nweights = 4*nn->q;
    make_vector(nn->weights, nn->nweights);
    for (int i = 0; i < nn->nweights; i++) {
        nn->weights[i] = ((double)rand())/RAND_MAX - 0.5;
    }
}


void Neural_Net_end(struct Neural_Net_PDE *nn) {
    free_vector(nn->weights);
}

void Neural_Net_eval(struct Neural_Net_PDE *nn, double x, double y) {
    int q     = nn->q;
    double *u = nn->weights;
    double *v = nn->weights +   q;
    double *w1 = nn->weights + 2*q;
    double *w2 = nn->weights + 3*q;

    for (int i = 0; i <= 2; i++) {
        for (int j = 0; j <= 2; j++) {
            nn->N[i][j] = 0.0;
        }
    }

    for (int i = 0; i < q; i++) {
        double z = u[i] + w1[i]*x+w2[i]*y;
        sigmoid(z, nn->sigma);
        // I could do double loops but i prefer this nice explicit form
        for (int ix = 0; ix <= 2; ix++) {
            for (int iy = 0; iy <= 2; iy++) {
                nn->N[ix][iy] += v[i]*pow(w1[i], ix)
                                     *pow(w2[i], iy)
                                     *nn->sigma[ix+iy];
            }
        }
        /*nn->N[0][0] += v[i]*nn->sigma[0];*/
        /*nn->N[1][0] += v[i]*pow(w1[i], 1)*nn->sigma[1];*/
        /*nn->N[2][0] += v[i]*pow(w1[i], 2)*nn->sigma[2];*/
        /*nn->N[0][1] += v[i]*pow(w2[i], 1)*nn->sigma[1];*/
        /*nn->N[0][2] += v[i]*pow(w2[i], 2)*nn->sigma[2];*/
        /*nn->N[1][1] += v[i]*pow(w1[i], 1)*/
                           /**pow(w2[i], 1)*nn->sigma[2];*/
    }
}

static double residual_at_x_y(struct Neural_Net_PDE *nn,
        double x, double y) {
    // Get N
    Neural_Net_eval(nn, x, y);
    // Get Phi
    nn->phi_func(nn, x, y);
    double u    = 
         nn->phi[0][0]*nn->N[0][0];
    double u_x  = 
         nn->phi[0][0]*nn->N[1][0]
        +nn->phi[1][0]*nn->N[0][0]; 
    double u_xx = 
         nn->phi[2][0]*nn->N[0][0]
        +nn->phi[0][0]*nn->N[2][0]
        +nn->phi[1][0]*nn->N[1][0]*2;
    double u_y  = 
         nn->phi[0][0]*nn->N[0][1]
        +nn->phi[0][1]*nn->N[0][0];
    double u_yy = 
         nn->phi[0][2]*nn->N[0][0]
        +nn->phi[0][0]*nn->N[0][2]
        +nn->phi[0][1]*nn->N[0][1]*2;
    double u_xy = 
         nn->phi[1][1]*nn->N[0][0]
        +nn->phi[0][0]*nn->N[1][1]
        +nn->phi[0][1]*nn->N[1][0]
        +nn->phi[1][0]*nn->N[0][1];
    return nn->PDE(x, y, u, u_x, u_y, u_xx, u_xy, u_yy);

}

double Neural_Net_residual(double *weights, int nweights, 
        void *params) {
    struct Neural_Net_PDE *nn = params;
    int nu = nn->nu; // the number of training points
    double sum = 0.0;

    for (int i = 0; i < nn->nweights; i++) {
        nn->weights[i] = weights[i];
    }

    for (int i =0; i < nu; i++) {
        double x = nn->training_points[i][0];
        double y = nn->training_points[i][1];
        double r = residual_at_x_y(nn, x, y);
        sum += r*r;
    }

    return sum;
}


double Neural_Net_error_vs_exact(struct Neural_Net_PDE *nn, 
        int m, int n) {
    if (nn->exact_sol == NULL) {
        fprintf(stderr,
                "unable to calculate the error since no exact solution is provided\n");
        exit(EXIT_FAILURE);
    }

    double max_err = 0.0;
    double a = nn->bb_xrange[0];
    double b = nn->bb_xrange[1];
    double c = nn->bb_yrange[0];
    double d = nn->bb_yrange[1];
    
    /*double x, y, z, err, sol;*/
    /*err = 0;*/
    double xh = (b-a)/m;
    double yh = (d-c)/n;
    for (int i = 0; i <= m; i++ ) {
        double x  = a + xh*i;
        for (int j = 0; j <= n; j++) {
            double y = c + yh*j;
            nn->phi_func(nn, x, y);
            if (nn->phi[0][0] < 0) {
                continue;
            }
            Neural_Net_eval(nn, x, y);
            double z = nn->N[0][0]*nn->phi[0][0];
            double sol = nn->exact_sol(x,y);
            double err = fabs(z - sol);
            if (err > max_err) {
                max_err = err;
            }
        }
    }
    return max_err;
}
