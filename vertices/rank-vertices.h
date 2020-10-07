#ifndef H_RANKVERTICES_H
#define H_RANKVERTICES_H

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

//static inline void rank_vertices(double *y, int m, int *ia, int *iy, int *iz);

#endif // H_RANKVERTICES_H
