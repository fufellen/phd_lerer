#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifndef SCOMPLEX_DEFINED
 // typedef MKL_Complex16 scomplex;
#define scomplex MKL_Complex16 
#endif

#include "mkl.h"
#include "mkl_blas.h"
#include "mkl_cblas.h"
#include "mkl_types.h"
#include "mkl_vml.h"


scomplex ** AllocC(int arow, int acol);
void DeAllocC(scomplex **p);
void ** AllocV(int arow, int acol, int size);
void DeAllocV(void **p);

int LineIt(scomplex **ca, int n, scomplex *b, double sigma);

int CLineEqC(scomplex **ca, int n, scomplex **cb, int nb);
int DLineEqD(double **ca, int n, double **cb, int nb);
int DLineEqFactor(int n, double **ca,int *ipiv);
int DLineEqSolve(int n,int nb,double **ca,double **cb,int * ipiv);
int CLineEqFactor(int n, scomplex **ca,int *ipiv);
int CLineEqSolve(int n,int nb,scomplex **ca,scomplex **cb,int * ipiv);	

//extern void    zgetrf(int *m,int *n, scomplex *a,int *lda,int *ipiv,int *info);
//extern void    zgetrs(char *trans,int *n,int *nrhs, scomplex *a,int *lda,int *ipiv, scomplex *b,int *ldb,int *info);
//extern void    dgetrf(int *m,int *n, double *a,int *lda,int *ipiv,int *info);
//extern void    dgetrs(char *trans,int *n,int *nrhs, double *a,int *lda,int *ipiv, double *b,int *ldb,int *info);
//extern void	   zcopy(int *n,scomplex *x,int *incx,scomplex *y,int *incy); 
//extern void    zgemv(char * trans, int *m, int *n, scomplex *alpha, scomplex *a, int *lda, scomplex *x, int *incx, scomplex *beta,scomplex *y, int *incy );


// операции с матрицами элементов double;

void MatrixSumD(int m,int n, double a, double **x, double **y);
void MatrixCopyD(int m,int n, double **x, double **y);
void MatrixDotD(int m,int n,  double **x, double **y, double *res);
void MatrixAbsD(int m,int n, double **x,  double *res);
void MatrixDivD(int n, double *x,  double *y,  double *res);
void MatrixScaleD(int m, int n, double **x,  double *y,  double **res);
void MatrixNormalizeD(int m, int n, double **x,  double *y,  double **res);
void Area(int n, double *x,double *y,double *s);
void VectorScaleD(int m, int n,  double a, double **x);
double VectorMaxAbs(int m, int n, scomplex **x);
void VectorScaleC(int m, int n,  scomplex a, scomplex **x);
double snorm(scomplex x);
void MatrixCopyC(int m,int n, scomplex **x, scomplex **y); //копируем x в y
int MatrixPack(scomplex ** A, // исходная матрица
			int n,			// размер матрицы
			double tol,		// точность определения нуля
			scomplex * values,// выходной массив со значениями ненулевых элементов матрицы
			int  * columns,		// выходной массив с номерами столбцов ненулевых элементов
			int * rowindex);	// выходной массив с номерами первых ненулевых элементов в массиве values

void ConjGrNR(scomplex **A,	// матрица системы (с точки зрения BLAS она транспонирована At)
			  int n,		// размер системной матрицы
			  scomplex ** B,// матрица правых частей (с точки зрения BLAS все нормально)
			  int m,		// число правых частей
			  scomplex ** G,// матрица начальных значений размер такой же как у B
			  double tol);	// точность контроля сходимости

void ConjGrNRJacobi(scomplex **A,	// матрица системы (с точки зрения BLAS она транспонирована At)
			  int n,		// размер системной матрицы
			  scomplex ** B,// матрица правых частей (с точки зрения BLAS все нормально)
			  int m,		// число правых частей
			  scomplex ** G,// матрица начальных значений размер такой же как у B
			  double tol);	// точность контроля сходимости
void ConjGrM(scomplex **A,	// матрица системы (с точки зрения BLAS она транспонирована At)
			  int n,		// размер системной матрицы
			  scomplex ** B,// матрица правых частей (с точки зрения BLAS все нормально)
			  int m,		// число правых частей
			  scomplex ** G,// матрица начальных значений размер такой же как у B
			  double tol);	// точность контроля сходимости
void BiCGStab(scomplex **A,	// матрица системы (с точки зрения BLAS она транспонирована At)
			  int n,		// размер системной матрицы
			  scomplex ** B,// матрица правых частей (с точки зрения BLAS все нормально)
			  int m,		// число правых частей
			  scomplex ** G,// матрица начальных значений размер такой же как у B
			  double tol);	// точность контроля сходимости
void DILUPivots(int n,			// размер матрицы системы
				scomplex ** A,	// матрица системы
				int w,			// число боковых диагоналей используемых для генерции вектора диагональных элементов
				scomplex * d);	// вектор содержащий обратные величины диагональных элементов DILU preconditioner

void DILUSolve(int n,			// размер матрицы системы
				scomplex ** A,	// матрица системы
				int w,			// число боковых диагоналей используемых для генерции вектора диагональных элементов
				scomplex * d,	// вектор содержащий обратные величины диагональных элементов DILU preconditioner
				scomplex * y,	// вектор правой части системы
			   scomplex *x);		// вектор решений
void PrintMatrixC (scomplex ** A, int n,int m, FILE * f);
void BiCGStabPreCond(scomplex **A,	// матрица системы (с точки зрения BLAS она транспонирована At)
			  int n,		// размер системной матрицы
			  scomplex ** B,// матрица правых частей (с точки зрения BLAS все нормально)
			  int m,		// число правых частей
			  scomplex ** G,// матрица начальных значений размер такой же как у B
			  double tol);	// точность контроля сходимости