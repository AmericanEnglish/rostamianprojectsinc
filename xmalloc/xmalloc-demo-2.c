#include <stdio.h>
#include "xmalloc.h"

int main(int argc, char *argv[]) {
    int GB = 1e9;
    int iter = 0;
    while (1) {
        void *x = xmalloc(GB);
        iter++;
        printf("%3d allocating %d bytes\n", iter, GB);
    }
    return 0;
}
