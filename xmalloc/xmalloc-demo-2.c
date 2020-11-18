#include <stdio.h>
#include "xmalloc.h"

int main(int argc, char *argv[]) {
    int GB = 1e9;
    int iter = 0;
    void *x;
    while (1) {
        x = xmalloc(GB);
        iter++;
        printf("%3d allocating %d bytes\n", iter, GB);
    }
    free(x);
    return 0;
}
