/* Замена того, что исходники VIBR_T брали из Intel MKL и библиотеки newlib+mat:
   LU-разложение и умножение матриц, векторные операции, ответ yesno() из стандартного ввода и libmat().
   Семантика zgetrf/zgetrs/zgemm - как у LAPACK/BLAS, хранение по столбцам. */
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "mkl_stub.h"

void Ag_(double lam, double *Re_eps, double *Im_eps);
void Au_(double lam, double *Re_eps, double *Im_eps);
void ZnO_(double lam, double *Re_eps, double *Im_eps);
void Si_(double lam, double *Re_eps, double *Im_eps);
void SiO2_(double lam, double *Re_eps, double *Im_eps);
void HfO_(double lam, double *Re_eps, double *Im_eps);
void AZO_(double lam, double *Re_eps, double *Im_eps);
void VO2_cold(double x, double *eps1, double *eps2);
void VO2_hot(double x, double *eps1, double *eps2);
double Cu(double x, int ri);

#ifdef FP_TRAP
#include <float.h>
__attribute__((constructor)) static void fp_trap(void)
{
	_controlfp(_controlfp(0, 0) & ~_EM_INVALID, _MCW_EM);
}
#endif

void zgetrf(int *m, int *n, MKL_Complex16 *a_, int *lda, int *ipiv, int *info)
{
	double complex *a = (double complex *)a_;
	int M = *m, N = *n, L = *lda, mn = M < N ? M : N;
	*info = 0;
	for (int k = 0; k < mn; k++) {
		int p = k;
		double big = cabs(a[k + k * L]);
		for (int i = k + 1; i < M; i++) {
			double v = cabs(a[i + k * L]);
			if (v > big) { big = v; p = i; }
		}
		ipiv[k] = p + 1;
		if (big == 0.0) {
			if (*info == 0) *info = k + 1;
			continue;
		}
		if (p != k)
			for (int j = 0; j < N; j++) {
				double complex t = a[k + j * L];
				a[k + j * L] = a[p + j * L];
				a[p + j * L] = t;
			}
		double complex d = 1.0 / a[k + k * L];
		for (int i = k + 1; i < M; i++) a[i + k * L] *= d;
		for (int j = k + 1; j < N; j++) {
			double complex f = a[k + j * L];
			if (f == 0.0) continue;
			double complex *cj = a + j * L, *ck = a + k * L;
			for (int i = k + 1; i < M; i++) cj[i] -= ck[i] * f;
		}
	}
}

void zgetrs(char *trans, int *n, int *nrhs, MKL_Complex16 *a_, int *lda, int *ipiv, MKL_Complex16 *b_, int *ldb, int *info)
{
	double complex *a = (double complex *)a_, *b = (double complex *)b_;
	int N = *n, L = *lda, LB = *ldb;
	char t = *trans;
	*info = 0;
	for (int r = 0; r < *nrhs; r++) {
		double complex *x = b + r * LB;
		if (t == 'N' || t == 'n') {
			for (int k = 0; k < N; k++) {
				int p = ipiv[k] - 1;
				if (p != k) { double complex s = x[k]; x[k] = x[p]; x[p] = s; }
			}
			for (int j = 0; j < N; j++) {
				double complex f = x[j];
				if (f == 0.0) continue;
				for (int i = j + 1; i < N; i++) x[i] -= a[i + j * L] * f;
			}
			for (int j = N - 1; j >= 0; j--) {
				x[j] /= a[j + j * L];
				double complex f = x[j];
				if (f == 0.0) continue;
				for (int i = 0; i < j; i++) x[i] -= a[i + j * L] * f;
			}
		} else {
			int cj = (t == 'C' || t == 'c');
			for (int j = 0; j < N; j++) {
				double complex s = x[j];
				for (int i = 0; i < j; i++) s -= (cj ? conj(a[i + j * L]) : a[i + j * L]) * x[i];
				x[j] = s / (cj ? conj(a[j + j * L]) : a[j + j * L]);
			}
			for (int j = N - 1; j >= 0; j--) {
				double complex s = x[j];
				for (int i = j + 1; i < N; i++) s -= (cj ? conj(a[i + j * L]) : a[i + j * L]) * x[i];
				x[j] = s;
			}
			for (int k = N - 1; k >= 0; k--) {
				int p = ipiv[k] - 1;
				if (p != k) { double complex s = x[k]; x[k] = x[p]; x[p] = s; }
			}
		}
	}
}

static double complex elem(const double complex *a, int L, char t, int i, int j)
{
	if (t == 'N' || t == 'n') return a[i + j * L];
	if (t == 'T' || t == 't') return a[j + i * L];
	return conj(a[j + i * L]);
}

void zgemm(char *ta, char *tb, int *m, int *n, int *k, MKL_Complex16 *al_, MKL_Complex16 *a_, int *lda,
           MKL_Complex16 *b_, int *ldb, MKL_Complex16 *be_, MKL_Complex16 *c_, int *ldc)
{
	double complex al = *(double complex *)al_, be = *(double complex *)be_;
	double complex *a = (double complex *)a_, *b = (double complex *)b_, *c = (double complex *)c_;
	int M = *m, N = *n, K = *k, LC = *ldc;
	for (int j = 0; j < N; j++) {
		for (int i = 0; i < M; i++) c[i + j * LC] = (be == 0.0) ? 0.0 : be * c[i + j * LC];
		for (int l = 0; l < K; l++) {
			double complex f = al * elem(b, *ldb, *tb, l, j);
			if (f == 0.0) continue;
			for (int i = 0; i < M; i++) c[i + j * LC] += elem(a, *lda, *ta, i, l) * f;
		}
	}
}

void vdDiv(int n, const double *a, const double *b, double *r) { for (int i = 0; i < n; i++) r[i] = a[i] / b[i]; }
void vdPowx(int n, const double *a, double b, double *r) { for (int i = 0; i < n; i++) r[i] = pow(a[i], b); }
void vdSqrt(int n, const double *a, double *r) { for (int i = 0; i < n; i++) r[i] = sqrt(a[i]); }

/* Эти функции нужны только итерационным решателям MLAP.C (ConjGr*, BiCGStab*), расчёт идёт через LU.
   Если вызов всё же произойдёт, программа остановится, а не даст неверное число. */
#define NOT_USED(name) void name(void) { printf("\n %s: не реализовано в port_c\n", #name); exit(1); }
NOT_USED(zcopy) NOT_USED(zaxpy) NOT_USED(zdotc) NOT_USED(dznrm2) NOT_USED(dnrm2)
NOT_USED(zscal) NOT_USED(zladiv) NOT_USED(zlacgv) NOT_USED(zgemv) NOT_USED(zgbmv)
NOT_USED(dgetrf) NOT_USED(dgetrs) NOT_USED(dcopy) NOT_USED(dscal) NOT_USED(DAXPY)

/* В исходной библиотеке ответ читается с клавиатуры, здесь из стандартного ввода, чтобы запускать из сценария. */
int yesno(void)
{
	char c = 'n';
	if (scanf(" %c", &c) != 1) return 0;
	return c == 'y' || c == 'Y';
}

/* Номера is - как в BIBL_mat.c и scripts/vibr_materials.py. Таблицы 10-19 в исходниках отсутствуют. */
void libmat(int is, double lam, double *Re_eps, double *Im_eps)
{
	switch (is) {
	case 1: *Re_eps = Cu(lam, 0); *Im_eps = Cu(lam, 1); break;
	case 2: Ag_(lam, Re_eps, Im_eps); break;
	case 3: Au_(lam, Re_eps, Im_eps); break;
	case 4: ZnO_(lam, Re_eps, Im_eps); break;
	case 5: Si_(lam, Re_eps, Im_eps); break;
	case 6: SiO2_(lam, Re_eps, Im_eps); break;
	case 7: HfO_(lam, Re_eps, Im_eps); break;
	case 8: AZO_(lam, Re_eps, Im_eps); break;
	case 20: VO2_cold(lam, Re_eps, Im_eps); break;
	case 21: VO2_hot(lam, Re_eps, Im_eps); break;
	default:
		printf("\n libmat: материала %d нет в исходниках\n", is);
		exit(1);
	}
}
