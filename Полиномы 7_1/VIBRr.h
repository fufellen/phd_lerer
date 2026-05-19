
#include <stdio.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <conio.h>
//#include <complex.h>

#include "newlib.h"
#include "matrix_4.h"
//#include "mlap.h"

#pragma warning (disable: 4996)
//#pragma warning (disable: LNK4229)

#define ESCape 0x1b

void VVOD(void);
void F_df_2(int sp, double *eps, double *mu, double *R, double *T);
void F_df_2_Re(int sp, double *eps, double *mu, double *R, double *T);
void F_df_IE_MK(int PP, double *R, double *T);
void J_c(double am, double alf, double l, double *rA, double *iA);
void J_s(double am, double alf, double l, double *rA, double *iA);

double  J(int n, double z);
double  N(int n, double z);
double  I_ex(int n, double z);

double  Pn(int n, double z);//Градштейн с. 1040 //Полиномы Лежандра

void U_(int n, int nu, int cn, int cnu, double gam2, double *rA, double *iA);
void V_(int PP, int n, int nu, double alf, double *rV, double *iV);

//void Q_mu(int cnu, double x, double *rA, double *iA);
double Q_mu(int mu, double l, int zn, double x);
void E_ext(int p, int mn_, double *rR, double *iR, double *rTA, double *iTA);

void F_GRIN_dop(int p, double alf, double *rR, double *iR, double *rTA, double *iTA);
//void F_df_dop(int PP, int m, int m_u, int n, int nu, double alf, double *rA, double *iA);
//void V_U_dop(int PP, int n, int nu, double alf, double *rA, double *iA);

void YSH(int nu, double gam, double *u, double*v);

void R_T(double *rR_, double *iR_);

void F_df_dop(int PP, int m, int m_u, int n, int nu, double alf, double *rA, double *iA);
void Y(double gam2, int nu, double  *rYnu, double *iYnu);
void R_T_dop(int PP, int m_u, int nu, double alf, double *rR_, double *iR_, double *rT_, double *iT_);

void compl_cos_sin(double x, double y, double *Rec, double *Imc, double *Res, double *Ims);
void compl_sqrt(double x, double y, double *Re, double *Im);
exp_compl(double x, double y, double *re, double *ie);//exp(iz)

double bess_dr_real(int n, double x, double eps, double nu, int i);//Функция Бесселя Jn+nu(x)/ x^nu
double gamma(double x);

void U(int n, int nu, int cn, int cnu, double bn, double bnu, double gam2, double *rA, double *iA);
void V(int n, int nu, int cn, int cnu, double bn, double bnu, double gam2, double *rA, double *iA);
void  Imod_K(int n, double z, double *I, double *K);
int graf1, graf2, M_X, M_Y, AM_Y[4], AM_X[4], AN_g[4], N_g, TEST, MS, SX, VAR,PRINT,PP;

double Pi, Pi_2, Pi05, f, TET, k, k2, Z0, a, h, G_R, AF[4], EPS[3], MU[3], kx, d, l, L, H, h2, eps[4], mu__[4], epsG[4],muG[4],
h2G;
double Rek, Imk, gam, r_tau[4], i_tau[4], r_mukap, i_mukap, r_eta[5], i_eta[5], **RX, **IX, **RA, **IA, rR_ext, iR_ext, rT_ext, iT_ext;

//scomplex **Alap,**Blap;


FILE *in, *in2, *out, *out1, *out2, *in1;
char inf[8], inf2[8], outf[8], out1f[8], out2f[8], in1f[8];
clock_t str,start;
void tim ( clock_t str );

