#ifndef H_NEURAL_NETS_PDE_H
#define H_NEURAL_NETS_PDE_H

struct Neural_Net_PDE {
    double (*PDE) (double x, double y, double u,
            double u_x,               double u_y,
            double u_xx, double u_xy, double u_yy);
    void (*phi_func)(struct Neural_Net_PDE *nn, double x, double y);
    double bb_xrange[2];
    double bb_yrange[2];

    int q;  // units per layer
    int nu; // number of training points
    double **training_points; 
    double (*exact_sol)(double x, double y);

    // no users allowed beyond this point
    int nweights;
    double *weights;
    double sigma[3];
    double phi[3][3];
    double N[3][3];
};

void Neural_Net_init(struct Neural_Net_PDE *nn);
void Neural_Net_end(struct Neural_Net_PDE *nn);
void Neural_Net_eval(struct Neural_Net_PDE *nn,
        double x, double y);
void Neural_Net_plot_with_maple(struct Neural_Net_PDE *nn,
        int n, int m, char *outfile);
void Neural_Net_plot_with_matlab(struct Neural_Net_PDE *nn,
        int n, int m, char *outfile);
double Neural_Net_residual(double *weights, int nweights, 
        void *params);
double Neural_Net_error_vs_exact(struct Neural_Net_PDE *nn, int m, 
        int n);


#endif // H_NEURAL_NETS_PDE_H
