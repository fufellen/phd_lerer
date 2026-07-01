
#include <stdio.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <conio.h>

#include "newlib.h"
#include "matrix_4.h"
//#include "mlap.h"

#pragma warning (disable: 4996)
//#pragma warning (disable: LNK4229)

#define ESCape 0x1b

void VVOD ( void );
void F_Compl( double *z,double *f);
void C_sqrt(double x,double y, double *u,double *v);
void C_exp(double x,double y, double *u,double *v);
void C_sh(double x,double y, double *u,double *v);
void C_ch(double x,double y, double *u,double *v);
void C_sh_ch(double x,double y, double *us,double *vs, double *uc,double *vc);
void Cn_sh_ch(double rg_,double ig_,double z,double h, double *us,double *vs, double *uc,double *vc);
void sq_compl(double x,double y,double *u,double *v);
void sq_compl_(double x,double y,double *u,double *v);
void SH_CTH(double x,double y,double *us,double *vs,double *uc,double *vc);
void sin_cos_C(double x,double y,double *us,double *vs,double *uc,double *vc);
void sin_C(double x,double y,double *us,double *vs);
void sinh_C(double x,double y,double *us,double *vs);
void cosh_C(double x,double y,double *us,double *vs);
void CONST(int N_d,int i,double Ro, double *te,double *qu,double *fi1,double *fi3,double *del,double *S1,double *S3);
double Max(double x, double y);
double Ag(double x, int ri);
double Au(double x, int ri);
double Cu(double x, int ri);
void Ag_(double lam, double *Re_eps, double *Im_eps);
void Au_(double lam, double *Re_eps, double *Im_eps);
void mult_a( double *x,double *y,double *z );
void div_a( double *x,double *y,double *z );
void General_COMPL ( int Nsl,double *X, int N_kor,double fto,double *bet_r,double *bet_i );
void F_1( void );

int graf2,graf1,Nsl,TV,ik,bum1,bum2,Pot,N_arg,Nkor,is[21][21],is_sub[5],E,jg,Nx,Ny,Nx_,nsub,Razb,n_cor,nc,hru;
	
double k,k2,Z0,AF[4],f,Pi_2,Pi,fto,Abet[6],tgd,EPS_MAX,EPS_min,EPS_r_[21],EPS_r[21],EPS_i[21],EPS_i_[21],
eps_r_[21][21],eps_r[21][21],eps_i[21][21],eps_i_[21][21],eps_sub_r[5],eps_sub_i[5],A_[21],B_[21],
	H_[21],DK[4][117],Ffu,X[7],k2_r[21],k2_i[21],k_r[21],SQ_05,pi2,n_[21],k_[21];

FILE *in,*out,*out1,*out2,*in1;
char inf[8], outf[8], out1f[8], out2f[8], in1f[8];
clock_t str,start;
void tim ( clock_t str );

