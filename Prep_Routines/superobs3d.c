/* File:          superobs3.c
 *
 * Created:       17 Nov 2009
 *
 * Last modified: 17 Nov 2009
 * Author:        Pavel Sakov
 *                NERSC
 *
 * Purpose:       Sorting of observations according to model grid cells.
 *
 * Description:   This file is an extension of the superobs.c for the 3D case.
 *
 * Modifications: none
 */

#include <math.h>
#include "cfortran.h"

#define IMAX 4096
#define JMAX 4096
#define KMAX 150

typedef struct {
    int i;
    int j;
    int k;
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
    else if (v1->k > v2->k)
	return 1;
    else if (v1->k < v2->k)
	return -1;
    return 0;
}

void sortgriddedobs3d(double pn, int ipiv[], int jpiv[], int kpiv[], int sorted[])
{
    int n = (int) pn;
    indexedvalue* iv = malloc(n * sizeof(indexedvalue));
    int i;

    for (i = 0; i < n; ++i) {
	int ii = ipiv[i];
	int jj = jpiv[i];
	int kk = kpiv[i];

	if (ii <= 0 || ii > IMAX || jj <= 0 || jj > JMAX || kk < 0 || kk > KMAX) {
	    fprintf(stderr, "ERROR: superobs.c: sortgriddedobs(): ipiv(%d) = %d or jpiv(%d) = %d or kpiv(%d) = %d out of bounds\n", i, ii, i, jj, i, kk);
	    exit(1);
	}
	iv[i].i = ii;
	iv[i].j = jj;
	iv[i].k = kk;
	iv[i].index = i;
    }
    
    qsort(iv, n, sizeof(indexedvalue), comp);

    for (i = 0; i < n; ++i)
	sorted[i] = iv[i].index + 1;

    free(iv);
}

FCALLSCSUB5(sortgriddedobs3d, SORTGRIDDEDOBS3D, sortgriddedobs3d, DOUBLE, PINT, PINT, PINT, PINT)
