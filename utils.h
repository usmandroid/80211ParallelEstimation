#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <complex.h>

#define PI M_PI
#define TWOPI (2.0*PI)
#define SIZESYMBOL 53*15
#define SAMPUTIL 53
#define OFDMBLK 15
#define P0 6
#define P1 20
#define P2 34
#define P3 48

void hermitian(long double complex **M, int row, int col, long double complex **res);
void multiply(long double complex **M1, int row1, int col1, long double complex **M2, int row2, int col2, long double complex **res);
void identity(long double complex **Identity, int size, double scalar);
void addition(long double complex **M1, int row1, int col1, long double complex **M2, int row2, int col2, long double complex **res);
void inverse(long double complex **A, int order, long double complex **Y);
int GetMinor(long double complex **src, long double complex **dest, int row, int col, int order);
long double complex CalcDeterminant( long double complex **mat, int order);
void fft_impl(double data[], int nn, int isign);
long double complex determinant_impl( long double complex **mat, int order);