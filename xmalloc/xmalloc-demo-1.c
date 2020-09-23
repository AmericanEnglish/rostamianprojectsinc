#include <stdio.h>
#include "xmalloc.h"

int main(int argc, char *argv[]) {
    int bytes[] = {1000, 0};
    for (int i = 0; i < 2; i++) {
        void *x = xmalloc(bytes[i]);
        printf("allocating %d bytes\n", bytes[i]);
        free(x);
        printf("memory freed\n");
    }
    return 0;
}
