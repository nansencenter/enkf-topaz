/* File:          superobs.c
 *
 * Created:       2 Sep 2008
 *
 * Last modified: 2 Sep 2008
 * Author:        Pavel Sakov
 *                NERSC
 *
 * Purpose:       Sorting of observations according to model grid cells.
 *
 * Description:   Given array of pivot indices for each observation, sort them
 *                in such a way that observations within each model grid cell
 *                will cluster together.
 *
 * Modifications: none
 */

#include <math.h>
#include "cfortran.h"

#define IMAX 4096
#define JMAX 4096

typedef struct {
    int i;
    int j;
    int index;
} indexedvalue;

static int comp(const void* p1, const void* p2)
{
    indexedvalue* v1 = (indexedvalue*) p1;
    indexedvalue* v2 = (indexedvalue*) p2;

    if (v1->i > v2->i)
	return 1;
    else if (v1->i < v2->i)
	return -1;
    else if (v1->j > v2->j)
	return 1;
    else if (v1->j < v2->j)
	return -1;
    return 0;
}

void sortgriddedobs(double pn, int ipiv[], int jpiv[], int sorted[])
{
    int n = (int) pn;
    indexedvalue* iv = malloc(n * sizeof(indexedvalue));
    int i;

    for (i = 0; i < n; ++i) {
	int ii = ipiv[i];
	int jj = jpiv[i];

	if (ii <= 0 || ii > IMAX || jj <= 0 || jj > JMAX) {
	    fprintf(stderr, "ERROR: superobs.c: sortgriddedobs(): ipiv(%d) = %d or jpiv(%d) = %d out of bounds\n", i, ii, i, jj);
	    exit(1);
	}
	iv[i].i = ii;
	iv[i].j = jj;
	iv[i].index = i;
    }
    
    qsort(iv, n, sizeof(indexedvalue), comp);

    for (i = 0; i < n; ++i)
	sorted[i] = iv[i].index + 1;

    free(iv);
}

FCALLSCSUB4(sortgriddedobs, SORTGRIDDEDOBS, sortgriddedobs, DOUBLE, PINT, PINT, PINT)
