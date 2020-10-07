#include "hilbert_matrix.h"
#include "array.h"

double **hilbert_matrix(int n) {
    // Get memory
    double **matrix;
    make_matrix(matrix, n, n);
    // Fill memory
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            matrix[i][j] = ((double) 1)/(i + j + 1);
        }
    }
    return matrix;
}

int main(int argc, char *argv[]) {
    int n = 8;
    double **H = hilbert_matrix(n);
    // Print?
    print_matrix("%7.3f", H, n, n);
    free_matrix(H);
    return 0;
}
