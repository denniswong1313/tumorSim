#ifndef ALLOC_H
#define ALLOC_H

double ** alloc2DContig(int nrows, int ncols);

int ** alloc2DContigInt(int nrows, int ncols);

void free2DContig(int nrows, int ncols, double ** array);

void free2DContigInt(int nrows, int ncols, int ** array);

#endif
