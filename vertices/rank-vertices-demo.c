/* rank-vertices-test.c
 * 
 * A program to test the rank_vertices() function (a part of nelder-mead).
 * 
 * It generates all vectors of length n with entries coming from
 * the set {0,1,2,...,n-1}.
 * 
 * Then it lets n vary from 2 to N, where N is #defined below.
 * 
 * For each generated vector, it runs the given rank_vertices() function.
 * 
 * If the result of rank_vertices() is incorrect, it prints a diagnostic,
 * otherwise no output is generated.
 * 
 * Rouben Rostamian
 * November 2003: n >= 3
 * November 2012: Extended to n = 2
 * March 2013: Added extra tests in test_result()
*/

#include <stdio.h>
#include "rank-vertices.h"
#define N 5

static void print_row(double *a, int n)
{
	int i;

	for (i = 0; i < n; i++)
		printf("%g ", a[i]);
}

static int incr(double *a, int n, int k)
{
	if (k == n)
		return 1;

	if (a[k] < n-1) {
		a[k]++;
		return 0;
	} else {
		a[k] = 0;
		return incr(a, n, k+1);
	}
}

/* Receives a vector a[] of length n containing entries from the set 0..(n-1).
   It "increments" the vector, as if the vector entries were the digits of
   an number in base n.  Returns 0 on success or 1 on failure, which happens
   when all vector entries have hit their max.
*/
static int increment(double *a, int n)
{
	return incr(a, n, 0);
}

static void test_result(double *a, int n, int ia, int iy, int iz)
{
	int i;

	/* are the indices distinct? */
	if ( (n == 2 && iy == iz)
			|| (n > 2 && (ia == iy || iy == iz || iz == ia))) {
		printf("vector: ");
		print_row(a, n);
		printf(" => ");
		printf("ia = %u, iy = %u, iz = %u\n", ia, iy, iz);
		return;
	}

	/* are the values ordered? */
	if (!(a[ia] <= a[iy] && a[iy] <= a[iz])) {
		printf("vector: ");
		print_row(a, n);
		printf("ia = %u, iy = %u, iz = %u\n", ia, iy, iz);
		printf("a[%u] <= a[%u] && a[%u] <= a[%u] failed!\n",
				ia, iy, iy, iz);
		return;
	}

	/* is ia a true min? */
	for (i = 0; i < n; i++)
		if (i != ia && a[i] < a[ia]) {
			printf("ia = %u but a[%u] < a[%u]!\n", ia, i, ia);
			return;
		}

	/* is iz a true max? */
	for (i = 0; i < n; i++)
		if (i != iz && a[i] > a[iz]) {
			printf("iz = %u but a[%u] > a[%u]!\n", iz, i, iz);
			return;
		}

	/* is iy correct? */
	for (i = 0; i < n; i++)
		if (i != iz && i != iy && a[i] > a[iy]) {
			printf("iy = %u but a[%u] > a[%u]!\n", iy, i, iy);
			return;
		}
}

int main(void)
{
	double a[N+1];
	int ia, iy, iz;
	int i, n;

	for (n = 2; n <= N; n++) {
		for (i = 0; i < n; i++)
			a[i] = 0;
		do {
			rank_vertices(a, n, &ia, &iy, &iz);
			test_result(a, n, ia, iy, iz);
		} while (increment(a, n) != 1);
	}

	return 0;
}
