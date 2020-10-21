#include <stdio.h>
#include "nelder-mead.h"
#include "math.h"

// Because why not? :D
#define eps 1.0e-3
#define XorE(X)((X < eps) ? eps : X)

static inline double tnorm2D(double x, double y) {
    return sqrt(pow(x, 2) + pow(y, 2));

}
static inline double W(double x) {
    return pow(x, 4)/24 +1/(12*pow(x,2)) - 1.0/8;
}
static double energy(double *x, int n, void *params) {
    // Assume params is just an array of doubles...
    double a = ((double *) params)[0];
    double b = ((double *) params)[1];
    return sqrt(2)*W(tnorm2D(x[0], x[1])/sqrt(2)) +
        W(tnorm2D(x[0]-1, x[1])) - a*(x[0]-1) - b*(x[1]-1);
}
/*static double obj_func(double*x, int n, void *params) {*/
    /*return 0;*/
/*}*/

int main(void){
    // Case 1
    /*double x[] = { 1.0, 1.0 };//the initial point */
    double x[] = { -1.0, 1.0 };//the initial point 
    double params[] = {0.0, -0.3};
    int evalcount;
    struct nelder_mead NM1 = {//Alert! C99-style initialization!
        .f      = energy,//the objective function
        .n      = 2, //the dimension of the space
        .s      = NULL, //delegate the construction of s
        .x      = x, //initial point/final point
        .h      = 1.0, //problem’s scale
        /*.h      = 0.1, //problem’s scale*/
        .tol    = 1.0e-3,//tolerance
        /*.tol    = 1.0e-4,//tolerance*/
        .maxevals = 1000,//cap on function evaluations
        .params = params,//no parameters
    };
    evalcount = nelder_mead(&NM1);
    printf("Case study 1\n");
    printf("Expected answer:   min = %13.10g at (%13.10g, %13.10g)\n",
            -0.6395281605, 0.8319642234, -1.2505278260);
    printf("Computed solution: min = %13.10g at (%13.10g, %13.10g)\n", NM1.minval, x[0], x[1]);
    printf("Converged after %d function evaluations\n", evalcount);

    double D1[]      = { 1.0, -1.0 };//the initial point 
    /*double D1[]      = { -1.0, 1.0 };//the initial point */
    double params2[] = { 0.0, -0.1};
    struct nelder_mead NM2 = {//Alert! C99-style initialization!
        .f      = energy,//the objective function
        .n      = 2, //the dimension of the space
        .s      = NULL, //delegate the construction of s
        .x      = D1, //initial point/final point
        .h      = 1.0, //problem’s scale
        /*.h      = 0.1, //problem’s scale*/
        .tol    = 1.0e-3,//tolerance
        /*.tol    = 1.0e-4,//tolerance*/
        .maxevals = 1000,//cap on function evaluations
        .params = params2,//no parameters
    };
    evalcount = nelder_mead(&NM2);
    printf("\n\nCase study 2, minimum at D1\n");
    printf("Expected answer:   min = %13.10g at (%13.10g, %14.11g)\n",
            -0.2047379473, 0.9208762303, -1.0926128103);
    printf("Computed solution: min = %13.10g at (%13.10g, %14.11g)\n", NM2.minval, D1[0], D1[1]);
    printf("Converged after %d function evaluations\n", evalcount);

    double D2[]      = { 1.0,  1.0};//the initial point 
    /*double D2[]      = { 1.0,  1.0};//the initial point */
    double params3[] = { 0.0, -0.1};
    struct nelder_mead NM3 = {//Alert! C99-style initialization!
        .f      = energy,//the objective function
        .n      = 2, //the dimension of the space
        .s      = NULL, //delegate the construction of s
        .x      = D2, //initial point/final point
        .h      = 1.0, //problem’s scale
        /*.h      = 0.1, //problem’s scale*/
        .tol    = 1.0e-3,//tolerance
        /*.tol    = 1.0e-4,//tolerance*/
        .maxevals = 1000,//cap on function evaluations
        .params = params3,//no parameters
    };
    evalcount = nelder_mead(&NM3);
    printf("\n\nCase study 2, minimum at D2\n");
    printf("Expected answer:   min = %13.10g at (%13.10g, %14.11g)\n",
            -0.0055676481, 1.1554672673, 0.8776603419);
    printf("Computed solution: min = %13.8g at (%13.10g, %14.11g)\n", NM3.minval, D2[0], D2[1]);
    printf("Converged after %d function evaluations\n", evalcount);
    return 0;
}
