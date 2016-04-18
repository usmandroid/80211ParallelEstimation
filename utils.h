#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <complex.h>
#include <mpi.h>

#define PI M_PI
#define TWOPI (2.0*PI)
#define SIZESYMBOL 53*15
#define SAMPUTIL 53
#define OFDMBLK 15
#define P0 5
#define P1 19
#define P2 33
#define P3 47

typedef struct {
    int numprocs, rank, res_LS;
    int tag1, tag2, tag3;
    long double complex conj1, conj2;
    MPI_Status status;
} Common_LT;

typedef struct {
    int numprocs, rank, res_LS;
    int tag1, tag2, tag3, tag4, tag5, tag6;
    long double complex H_PILOTS[4];
    long double H_PILOTS_real[4], H_PILOTS_imag[4];
    MPI_Status status;
} Common_PS;

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

int GetMinor(long double complex **src, long double complex **dest, int row, int col, int order);
long double complex CalcDeterminant( long double complex **mat, int order);
long double complex determinant_impl( long double complex **mat, int order);
long double complex determinant_impl_rec( long double complex **mat, int order);
long double complex determinant_impl_omp( long double complex **mat, int order);
void swap_rows(long double complex **mat, int order, int row1, int row2);
void swap_cols(long double complex **mat, int order, int row1, int row2);

void inverse(long double complex **A, int order, long double complex **Y);
void inverse_omp(long double complex **A, int order, long double complex **Y);
void inverse_mpi(long double complex **A, int order, long double complex **Y, Common_PS *commonPS, int argc, char *argv[]);
void inverse_mpi_beta(long double complex **A, int order, long double complex **Y, Common_PS *commonPS, int argc, char *argv[]);

void fft_impl(double data[], int nn, int isign);
double sinc(double input);

void printVect(long double complex *vec, int size,char *name);
void printMatrix(long double complex **mat, int rows, int cols, char *name);
void printMatrixReal(long double **mat, int rows, int cols, char *name);
void printMatrix2(long double complex **mat, int rows, int cols, char *name);

void test_master(int numprocs);
void test_slave(int rank);
void multiply_mpi(long double complex **M1, int row1, int col1, long double complex *vec, int col2, long double complex **res, int from, int to);

void complexToDouble(int order, long double complex **matrix, long double **matrix_Re, long double **matrix_Im);
void doubleToComplex(int order, long double complex **matrix, long double **matrix_Re, long double **matrix_Im);

int malloc2dLongDouble(long double ***array, int n, int m);
int malloc2dLongDoubleComplex(long double complex ***array, int n, int m);
int free2dLongDouble(long double ***array);
int free2dLongDoubleComplex(long double complex ***array);