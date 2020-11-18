#include "neural-nets-pde.h"
#include <math.h>
#include <stdio.h>

static double exact_sol(double x, double y) {
    return 16*x*(1-x)*y*(1-y);
}

static double my_pde(double x, double y, double u_x, 
        double u_y, double u_xx, double u_xy, double u_yy) {
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

}
