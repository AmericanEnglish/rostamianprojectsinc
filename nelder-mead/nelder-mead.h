#ifndef H_NELDERMEAD_H
#define H_NELDERMEAD_H

struct nelder_mead {
    double (*f)(double *x, int n, void *params);
    int n;
    double **s;
    double *x;
    double h;
    double tol;
    int maxevals;
    double minval;
    void *params;
};
//static inline void rank_vertices(double *y, int m, int *ia, int *iy, int *iz);
//static inline void rank_vertices(double *y, int m, int *ia, int *iy, int *iz);
//static void get_centroid(double **s, int n, int iz, double *C);
//static inline void transform(double *P, double *Q, int n, double beta, double *R);
//static void shrink(double **s, int n, int ia);
//static inline void replace_row(double **s, int i, double **row);
//static int done(double **s, int n, double *y, int ia, int iz, double err2);
int nelder_mead(struct nelder_mead *nm);
#endif // H_NELDERMEAD_H
