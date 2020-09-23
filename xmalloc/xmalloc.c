/*#include <stdio.h>*/
/*#include <stdlib.h>*/
#include "xmalloc.h"

void *malloc_of_exit(size_t nbytes, const char *file, int line) {
    void *x;
    if ((x = malloc(nbytes)) == NULL) {
        fprintf(stderr, "%s:line %d: malloc() of %zu bytes failed!\n",
                file, line, nbytes);
        exit(EXIT_FAILURE);
    }
    else {
        return x;
    }
}
