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

void hermitian(double complex **M,double complex **res);
void multiply(double complex **M1,double complex **M2,double complex **res);
void identity(double complex **Identity,int size,double scalar);
void addition(double complex **M1,double complex **M2,double complex **res);
void inverse(double complex **A, int order, double complex **Y);
int GetMinor(double complex **src, double complex **dest, int row, int col, int order);
double complex CalcDeterminant( double complex **mat, int order);
void fft_impl(double data[], int nn, int isign);
long double complex determinant_impl( double complex **mat, int order);