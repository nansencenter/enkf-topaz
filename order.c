/* File:          order.c
 *
 * Created:       2 Mar 2008
 *
 * Last modified: 2 Mar 2008
 * Author:        Pavel Sakov
 *                NERSC
 *
 * Purpose:       Put indices of an array of double in an order of increasing
 *                value.
 *
 * Description:   Given a double array x[n], sort its subset specified by an
 *                integer array of indices good[ngood] and return the indices
 *                of sorted elements in the integer array inorder[ngood].
 *                
 *                It is assumed that good[ngood] stores the "fortran" indices
 *                (from 1 to N rather than from 0 to N - 1).
 *
 * Modifications: none
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cfortran.h"

typedef struct {
    int index;
    double v;
} indexedvalue;

static int comp(const void* p1, const void* p2)
{
    indexedvalue* v1 = (indexedvalue*) p1;
    indexedvalue* v2 = (indexedvalue*) p2;

    if (v1->v > v2->v)
	return 1;
    else if (v1->v < v2->v)
	return -1;
    return 0;
}

/** Sorts a specified subset within an array of double according to values.
 *
 * Given a double array x[n], sorts its subset specified by an integer array
 * good[ngood] and returns the indices of sorted elements in the preallocated
 * integer array inorder[ngood].
 *
 * It is assumed that good[ngood] stores the "fortran" indices (from 1 to N
 * rather than from 0 to N - 1).
 *
 * @param pn Number of elements in the data array
 * @param x Data array
 * @param pngood Number of elements in the data array to be sorted
 * @param good Indices of the elements in the data array to be sorted
 * @param inorder Output array of size of `ngood' such that the corresponding
 *                elements of the data array are in increasing order
 */
void order(double pn[], double x[], double pngood[], int good[], int inorder[])
{
    int n = (int) pn[0];
    int ngood = (int) pngood[0];
    indexedvalue* iv = NULL;
    int i;
    
    if (n <= 0) {
	for (i = 0; i < ngood; ++i)
	    inorder[i] = -1;
	return;
    }

    iv = malloc(n * sizeof(indexedvalue));
    if (n < ngood) {
	fprintf(stderr, "ERROR: order.c: order(): size of the data = %d is less than the requested size of the sorted array %d\n", n, ngood);
	exit(1);
    }

    /*
     * a bit of quality control
     */
    for (i = 0; i < ngood; ++i) {
	double xx;

	if (good[i] < 1 || good[i] > n) {
	    fprintf(stderr, "ERROR: order.c: order(): good[%d] = %d, n = %d\n", i, good[i], n);
	    exit(1);
	}
	xx = x[good[i] - 1];
	if (isnan(xx) || fabs(xx) > 1.0e+10 || xx == -999.0) {
	    fprintf(stderr, "ERROR: order.c: order(): x[%d] = %.15g\n", good[i] - 1, xx);
	    exit(1);
	}
    }

    for (i = 0; i < ngood; ++i) {
	iv[i].index = good[i];
	iv[i].v = x[good[i] - 1];
    }

    qsort(iv, ngood, sizeof(indexedvalue), comp);

    for (i = 0; i < ngood; ++i)
	inorder[i] = iv[i].index;

    free(iv);
}

FCALLSCSUB5(order, ORDER, order, PDOUBLE, PDOUBLE, PDOUBLE, PINT, PINT)
