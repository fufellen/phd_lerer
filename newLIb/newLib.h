
#pragma once
#include "math.h"
#include "time.h"
#include "stdio.h"
#include "stdlib.h"
#pragma warning (disable: 4996)

/* Mc'Donald function of zero order          double KKK0(double x); */
/* Mc'Donald function of 1-th order           double KKK1(double x);*/
double K0( double x );
double K1( double x );
double Knp( int N,double x );/* Производная от Функции Макдональда */
double Kn( int N,double x );/* Функция Макдональда */

void Hermit( int N,double x,double *He );
void GEGEN( int N,double x,double lam,double *Ge );
double Arth( double x );
double th( double x );

/* modified Bessel function of zero order    double I0(double x);*/ 
/* modified Bessel function of first order    double I1(double x);*/
/* Mc'Donald function of zero order          double KKK0(double x); */
/* Mc'Donald function of 1-th order           double KKK1(double x);*/

/*
   Квадратурные узлы и коэффициенты квадратуры наивысшей точности для
   интералов на [0,1] от x*f(x). Крылов стр. 163
*/

void UZLI( int n,double *x,double *A );

/*
   Квадратурные узлы и коэффициенты квадратуры наивысшей точности Гаусса для
   интералов на [-1,1]. Крылов стр. */

void UZLI_GAU( int n,double *x,double *A );
void KVADR_GAU(int m, double a, double b, double *xk,  double *Ak );

double K(double x); /* эллиптический интеграл 1-го рода */

double I0(double x);
double I1(double x);
double KKK0(double x);
double KKK1(double x);
double Im( int m, double z );
double Imp( int m, double z );

void BESSEL(double x,double y, int n, double ftol,
			double *J, double *N);
			/* Функции Бесселя
					   1-го рода ( J[0] + i * J[1] ) и
					   2-го рода ( N[0] + i * N[1] )
					   порядка n
			   комплексного аргумента x + i y.
			   Разложение в ряд
			*/

int yesno();
int answer( char *code );
void mbess( int i, int ns, int n, double x[], double nu, double **j );
double HORKUZ( double a, double b, double fa, double fb,
			double (*funk)(double), double* ff, double tol);
double detsig( int n, double **d );
void CELINEG(int M, int N, double** RA, double** IA);

void detcompl( int n, double ** d , double ** id );
void mult( double x1,double y1,double x2,double y2,double *u,double *v );
void divide(double x,double y,double x1,double y1,double *u,double *v);

void *amat_1d(  /* Allocate one-dimensional array */
int lo,         /* Lower limit, inclusive */
int hi,         /* Upper limit, inclusive */
int el_size );   /* Elements size */
void fmat_1d(   /* Free one-dimensional array */
void *ptr,      /* Array pointer, allocated by amat_1d */
int lo,         /* Lower limit, inclusive */
int el_size )   /* Elements size */;
void **amat_2d( /* Allocate two-dimensional array */
int lo_row,     /* First lower limit, inclusive */
int hi_row,     /* First upper limit, inclusive */
int lo_col,     /* Second lower limit, inclusive */
int hi_col,     /* Second upper limit, inclusive */
int el_size )   /* Elements size */;
void fmat_2d(           /* Free two-dimensional array */
register void **ptr,    /* Array pointer, allocated by amat_2d */
int lo_row,             /* First lower limit, inclusive */
int hi_row,             /* First upper limit, inclusive */
int lo_col,             /* Second lower limit, inclusive */
int el_size )           /* Elements size */;
void *MemCalloc( size_t num, size_t size );         /* calloc */
void MemFree( void *memblock );
void *amem_1d(          /* Allocate one-dimensional array */
int dim,                /* Upper limit, exclusive */
int el_size );          /* Elements size */

void fmem_1d(           /* Free one-dimensional array */
void *ptr );            /* Array pointer, allocated by amem_1d */

void tim ( clock_t start );
void fSkipComment( FILE *stream );

void mult_a( double *x,double *y,double *z );//компл умножение
void div_a( double *x,double *y,double *z );//компл деление
void subt_a( double *x,double *y,double *z );//компл вычитание
void addit_a( double *x,double *y,double *z );//компл сложение
double mod( double *x,double *y);//модуль разности двух компл чисел
void NEW_C( double *a, double st, void func(double *z,double *f), double* ff, double tol);
void CEL( int n, int nm, double **a);
int gauss( int m, double* p, double* k );
int jacob1( int m, double* p, double* k );
double In( int n, double x, int N);

void ***amat_3d(    /* Allocate three-dimensional array */
int lo_depth,       /* First lower limit, inclusive */
int hi_depth,       /* First upper limit, inclusive */
int lo_row,         /* Second lower limit, inclusive */
int hi_row,         /* Second upper limit, inclusive */
int lo_col,         /* Third lower limit, inclusive */
int hi_col,         /* Third upper limit, inclusive */
int el_size )       /* Elements size */;
void fmat_3d(           /* Free three-dimensional array */
register void ***ptr,   /* Array pointer, allocated by amat_3d */
int lo_depth,           /* First lower limit, inclusive */
int hi_depth,           /* First upper limit, inclusive */
int lo_row,             /* Second lower limit, inclusive */
int hi_row,             /* Second upper limit, inclusive */
int lo_col,             /* Third lower limit, inclusive */
int el_size )           /* Elements size */;

void Sort( int m, double *x );
void SortInt( int m, int *x );
double GENERAL_L(double COM,double ST,double FN, double ftol,double (*F_df)(double),int *ik);
void NEW_KV_GAUS(int M, double a, double b, double *xk,double *Ak );
void cZeydel( int N, int It, double **RA, double **IA, double *RZ, double *IZ, double eps );
void ITERSLAU(int M, int It, double ftol, double **RA, double **IA,double ReC, double ImC, double *V );
void Legendr_P( int N,  double x, double **Pnm );
void Legendr( int N, double x, double *Pnm);