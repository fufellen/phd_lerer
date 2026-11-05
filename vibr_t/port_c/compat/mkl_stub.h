#pragma once
/* Объявления того немногого из Intel MKL, чем пользуются исходники VIBR_T. Реализация - vibr_compat.c. */
typedef struct { double real, imag; } MKL_Complex16;
void zgetrf(int *m, int *n, MKL_Complex16 *a, int *lda, int *ipiv, int *info);
void zgetrs(char *trans, int *n, int *nrhs, MKL_Complex16 *a, int *lda, int *ipiv, MKL_Complex16 *b, int *ldb, int *info);
void zgemm(char *ta, char *tb, int *m, int *n, int *k, MKL_Complex16 *alpha, MKL_Complex16 *a, int *lda,
           MKL_Complex16 *b, int *ldb, MKL_Complex16 *beta, MKL_Complex16 *c, int *ldc);
void vdDiv(int n, const double *a, const double *b, double *r);
void vdPowx(int n, const double *a, double b, double *r);
void vdSqrt(int n, const double *a, double *r);
