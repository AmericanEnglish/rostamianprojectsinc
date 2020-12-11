#ifndef H_PME_SOLVE_H
#define H_PME_SOLVE_H

struct pme_solve {
    double a;
    double b;
    double T;
    int n;
    int m;
    double (*ic)  (double x);
    double (*bcL) (double t);
    double (*bcR) (double t);
    double **u;
    double (*exact_sol)(double x, double t);
    double error;
    char *maple_out;
    char *matlab_out;
    char *geomview_out;
};

void show_usage_and_exit(char *progname);
void pme_solve(struct pme_solve *prob);

#endif // H_PME_SOLVE_H
