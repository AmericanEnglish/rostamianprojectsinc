#ifndef H_XMALLOC_H
#define H_XMALLOC_H

#include <stdio.h>
#include <stdlib.h>


void *malloc_of_exit(size_t nbytes, const char *file, int line);
// Easy automatic usage of the function
#define xmalloc(nbytes) malloc_or_exit((nbytes), __FILE__, __LINE__)

#endif
