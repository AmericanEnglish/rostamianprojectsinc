#ifndef H_NEURAL_NET_ODE_H
#define H_NEURAL_NET_ODE_H

struct Neural_Net_ODE {
    double (*ODE)(double x, double u, double u_x, double u_xx);
    double a;                      // left endpoint of the interval
    double b;                      // right endpoint of the interval
    int q;                         // number of units in the hidden layer
    int nu;                        // the number of training points
    double *training_points;       // the array of training points
    double (*exact_sol)(double x); // NULL if no exact solution avail

    // no user modifiable parts beyond this point
    int nweights;    // 3 x q
    double *weights; // the array of u, v, w
    double sigma[3]; // array to hold sig, dsig, ddsig
    double phi[3];   // array phi, dphi, ddphi
    double N[3];     // array to hold N, dN, ddN
}

void Neural_Net_init(struct Neural_Net_ODE *nn);
void Neural_Net_end(struct Neural_Net_ODE *nn);
void Neural_Net_eval(struct Neural_Net_ODE *nn);
void Neural_Net_phi(struct Neural_Net_ODE *nn);
void Neural_Net_plot_with_maple(struct Neural_Net_ODE *nn, int n, char *outfile);
void Neural_Net_plot_with_matlab(struct Neural_Net_ODE *nn, int n, char *outfile);
double Neural_Net_residual(double *weights, int nweights,
        void *params);
double Neural_Net_error_vs_exact(struct Neural_Net_ODE *nn, int n);

#endif /* H_NEURAL_NET_ODE_H */
