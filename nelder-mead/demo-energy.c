#include <stdio.h>
#include "nelder-mead.h"
#include "math.h"

// Because why not? :D
#define eps 1.0e-3
#define XorE(X)((X < eps) ? eps : X)
#define W(X) (pow(X, 4)/24.0 +1/(12.0*pow(X,2)) -1/8.0)
/*#define W(X) (pow(XorE(X), 4)/24.0 +1/(12.0*pow(XorE(X),2)) -1/8.0)*/
#define NORM2(X, Y) (pow(X, 2)+pow(Y, 2))

static double energy(double *x, int n, void *params) {
    // Assume params is just an array of doubles...
    double a = ((double *) params)[0];
    double b = ((double *) params)[1];
    return sqrt(2)*W(NORM2(x[0], x[1])/sqrt(2)) +
        W(NORM2(x[0]-1, x[1])) - a*(x[0]-1) - b*(x[1]-1);
}
/*static double obj_func(double*x, int n, void *params) {*/
    /*return 0;*/
/*}*/

int main(void){
    double x[] = { 1.0, 1.0 };//the initial point 
    double params[] = {0.0, -0.3};
    int evalcount;
    struct nelder_mead NM1 = {//Alert! C99-style initialization!
        .f      = energy,//the objective function
        .n      = 2, //the dimension of the space
        .s      = NULL, //delegate the construction of s
        .x      = x, //initial point/final point
        .h      = 0.1, //problem’s scale
        .tol    = 1.0e-4,//tolerance
        .maxevals = 1000,//cap on function evaluations
        .params = params,//no parameters
    };
    evalcount = nelder_mead(&NM1);
    if(evalcount > NM1.maxevals) {
        printf("No convergence after %d function evaluation\n",evalcount);
    }
    else{
        printf("converged after %d function evaluations\n",evalcount);
        printf("Computed solution: min = %g at (%g, %g)\n", NM1.minval, x[0], x[1]);
    }
    return 0;
}
