#include <stdio.h>
#include <stdlib.h>
#include <math.h>  

double ** alloc2DContig(int nrows, int ncols)
{
	int i, j;
	double **array;

	array = (double**)malloc(nrows * sizeof(double*));
	array[0] = (double*)calloc(sizeof(double), nrows*ncols);
	for (i = 1; i<nrows; i++)
	{
		array[i] = array[0] + i*ncols;
	}

	return array;
}

int ** alloc2DContigInt(int nrows, int ncols)
{
	int i, j;
	int **array;

	array = (int**)malloc(nrows * sizeof(int*));
	array[0] = (int*)calloc(sizeof(int), nrows*ncols);
	for (i = 1; i<nrows; i++)
	{
		array[i] = array[0] + i*ncols;
	}

	return array;
}

void free2DContig(int nrows, int ncols, double ** array)
{
	int i, j;

	for (i = 1; i<nrows; i++)
	{
		array[i] = NULL;
	}
	free(array[0]);
	free(array);
}

void free2DContigInt(int nrows, int ncols, int ** array)
{
	int i, j;

	for (i = 1; i < nrows; i++)
	{
		array[i] = NULL;
	}
	free(array[0]);
	free(array);
}