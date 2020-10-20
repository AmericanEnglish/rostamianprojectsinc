#include <stdio.h>

#include "array.h"
#include "nelder-mead.h"

#define REFLECT  1.0
#define EXPAND   2.0
#define CONTRACT 0.5 
#define SHRINK   0.5

//inline void rank_vertices(double *y, int m, int *ia, int *iy, int *iz)
static inline void rank_vertices(double *y, int m, int *ia, int *iy, int *iz) {
    // m >= 3
    // ia != iy != iz
    // m = 2
    // ia = iy != iz
    // ia = min y
    // iz = max y
    // iy = max y\{y[iy]} 
    
    // For easy typing, use these standins
    // last largest, second largest, first largest...
    int last, second, first;

    // Process the first two to satisfy the edge case
    if (y[0] <= y[1]) {
        last   = 0;
        second = 0;
        first  = 1;
    }
    else {
        last   = 1;
        second = 1;
        first  = 0;
    }
    // If m >= 3
    for (int i = 2; i < m; i++) {
        // If the first largest moves then the second largest moves too
        if (y[i] > y[first]) {
            second = first;
            first  = i;
        } 
        // In case the first index is the largest but not the second largest
        else if (i != first && i != last && y[i] >= y[second]) {
            second = i;
        }
        else if (y[i] < y[last]) {
            last   = i;
        }
    }
    // "return" values
    *ia = last;
    *iy = second;
    *iz = first;
}

static inline void get_centroid(double **s, int n, int iz, double *C) {
    // Note that the matrix is composed of row vectors
    int col, row;
    for (int col = 0; col < n ; col++) {
        C[col] = 0;
    }
    // We are skipping crow == iz, do the rows before iz
    for (row = 0; row < iz; row++) {
        for (col = 0; col < n; col++) {
            C[col] += s[row][col];
        }
    }
    for (row = (iz + 1); row < n+1; row++) {
        for (col = 0; col < n; col++) {
            C[col] += s[row][col];
        }
    }
    // Doing this after helps prevent round-off
    for (int col = 0; col < n ; col++) {
        C[col] = C[col] / ((double) n);
    }
}

static inline void transform(double *P, double *Q, int n, double beta, double *R) {
    int i;
    double oneMB = 1 - beta;
    for (i = 0; i < n; i++) {
        R[i] = oneMB*P[i] + beta*Q[i];
    }
}
static inline void shrink(double **s, int n, int ia) {
    int row;
    // P is the base point
    // R is the final point
    // Q is the current point
    // Shrink PQ to PR
    // => Shrink all vectors to ia
    for (row = 0; row < ia; row++) {
        transform(s[ia], s[row], n, SHRINK, s[row]);
    }
    for (row = (ia + 1); row < n+1; row++) {
        transform(s[ia], s[row], n, SHRINK, s[row]);
    }
}

static inline void replace_row(double **s, int i, double **row) {
    double *tmp = s[i];
    s[i] = *row;
    *row = tmp;
}

static inline int done(double **s, int n, double *y, int ia, int iz, double err2) {
    // ||s[iz] - s[ia]||^2 < eps^2
    // |y[iz] - y[ia]| <= eps^2
    // eps^2 == err2
    int col;
    double norm2 = 0;
    double diff;
    for (col = 0; col < n; col++) {
        diff   = s[iz][col] - s[ia][col];
        norm2 += diff*diff;
    }
    diff = y[iz] - y[ia];
    // To avoid using having to link the math library I plan to leverage the
    // following:
    // For strictly positive numbers if x < y => x^2 < y^2, so
    // ||s[iz]-s[ia]||^2 <= err2
    // And also convert from absolute value to a full inequality
    // |y[iz]-y[ia]| <= err2 => -err2 <= y[iz]-y[ia] <= err2
    return (norm2 < err2) && ((-err2 <= diff) && (diff <= err2));
}

int nelder_mead(struct nelder_mead *nm) {
    // 1st half
    double **s  = nm->s;
    int n       = nm->n;
    double h    = nm->h;
    double tol  = nm->tol;
    double err2 = (h*tol)*(h*tol);
    double *y, *C, *Pr, *Pe, *Pc;
    double yr, ye, yc;
    int ia, iy, iz;
    int simplex_to_be_freed = 0;
    int fevalcount;
    int row, col;

    make_vector(y, n+1);
    make_vector(Pr, n);
    make_vector(Pe, n);
    make_vector(Pc, n);
    make_vector(C, n);
    // 2nd half
    if (s == NULL) {
        make_matrix(s, n+1, n);
        simplex_to_be_freed = 1;
        // set values of s
        for (row = 0; row < n+1; row++) {
            for (col = 0; col < n; col++) {
                s[row][col] = nm->x[col];
                if ((row - 1) == col) {
                    s[row][col] += h;
                }
            }
        }
        /*print_matrix("%f ", s, n+1, n);*/
    }
    for (row = 0; row < n+1; row++) {
        y[row] = nm->f(s[row], n, nm->params);
    }
    fevalcount = n+1;
    while (fevalcount <= nm->maxevals) {
        rank_vertices(y, n+1, &ia, &iy, &iz);
        if (done(s, n, y, ia, iz, err2)) {
            nm->minval = y[ia];
            // copy best vertex into nm->x
            for (col = 0; col < n; col++) {
                nm->x[col] = s[ia][col];
            }
            break;
        }
        get_centroid(s, n, iz, C);

        /*printf("C = "); print_vector("%24.16f ", C, n);*/
        /*printf("s = \n"); print_matrix("%24.16f ", s, n+1, n);*/

        transform(C, s[iz], n, -REFLECT, Pr);
        yr = nm->f(Pr, n, nm->params);
        fevalcount++;

        if (yr < y[ia]) { // case 1
            transform(C, Pr, n, EXPAND, Pe);
            ye = nm->f(Pe, n, nm->params);
            fevalcount++;
            if (ye < yr) {
                replace_row(s, iz, &Pe);
                y[iz] = ye;
            } 
            else {
                replace_row(s, iz, &Pr);
                y[iz] = yr;
            }
        }
        else if (yr < y[iy]) { // case 2
            replace_row(s, iz, &Pr);
            y[iz] = yr; 
        }
        else { 
            if (yr < y[iz]) { //case 3
                //  R = (1-b)P+bQ
                transform(C, Pr, n, CONTRACT, Pc);
            // Could be an else, but for completion
            }
            else if (y[iz] <= yr) { // case 4
                transform(C, Pr, n, -CONTRACT, Pc);
            }
            yc = nm->f(Pc, n, nm->params);
            fevalcount++;
            if (yc < yr) {
                replace_row(s, iz, &Pc);
                y[iz] = yc;
            }
            else { // else, shrink
                shrink(s, n, ia);
                for (row = 0; row < ia; row++) {
                    y[row] = nm->f(s[row], n, nm->params);
                }
                for (row = ia; row < n+1; row++) {
                    y[row] = nm->f(s[row], n, nm->params);
                }
                fevalcount += n;
            }
        }
    }
    /*print_matrix("%f ", s, n+1, n);*/
    // Free all vectors
    free_vector(y);
    free_vector(Pr);
    free_vector(Pe);
    free_vector(Pc);
    free_vector(C);
    if (simplex_to_be_freed) {
        free_matrix(s);
    }
    return fevalcount;
}
