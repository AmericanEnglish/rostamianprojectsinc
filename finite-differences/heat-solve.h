#ifndef H_HEAT_SOLVE_H
#define H_HEAT_SOLVE_H

enum method {FD_undefine, FD_explicit, FD_implicit,
    FD_crank_nicolson, FD_seidman_sweep};

struct heat_solve {
    double a;
    double b;
    double T;
    int n;
    int m;
    double (*ic)  (double x);
    double (*bcL) (double t);
    double (*bcR) (double t);
    enum method method;
    double **u;
    double (*exact_sol) (double x, double t);
    double error;
    char *maple_out;
    char *matlab_out;
    char *geomview_out;
};

void show_usage_and_exit(char *progname);
void heat_solve(struct heat_solve *prob);

#endif
