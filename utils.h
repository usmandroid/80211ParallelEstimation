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
#define P0 5
#define P1 19
#define P2 33
#define P3 47

void hermitian(long double complex **M, int row, int col, long double complex **res);
void hermitian_omp(long double complex **M, int row, int col, long double complex **res);
void multiply(long double complex **M1, int row1, int col1, long double complex **M2, int row2, int col2, long double complex **res);
void multiply_omp(long double complex **M1, int row1, int col1, long double complex **M2, int row2, int col2, long double complex **res);
void multiplyVxVeqM(long double complex **M1, int row1, int col1, long double complex **M2, int row2, int col2, long double complex **res);
void multiplyVxVeqM_omp(long double complex **M1, int row1, int col1, long double complex **M2, int row2, int col2, long double complex **res);
void identity(long double complex **Identity, int size, double scalar);
void identity_omp(long double complex **Identity, int size, double scalar);
void addition(long double complex **M1, int row1, int col1, long double complex **M2, int row2, int col2, long double complex **res);
void addition_omp(long double complex **M1, int row1, int col1, long double complex **M2, int row2, int col2, long double complex **res);
void inverse(long double complex **A, int order, long double complex **Y);
void inverse_omp(long double complex **A, int order, long double complex **Y);
int GetMinor(long double complex **src, long double complex **dest, int row, int col, int order);
long double complex CalcDeterminant( long double complex **mat, int order);
void fft_impl(double data[], int nn, int isign);
long double complex determinant_impl( long double complex **mat, int order);
void printVect(long double complex *vec, int size,char *name);
void printMatrix(long double complex **mat, int rows, int cols, char *name);
double sinc(double input);