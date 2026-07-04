
#include <stdio.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <conio.h>

#include "newlib+mat.h"
#include "matrix_4.h"
#include "mlap.h"

#pragma once 

#pragma warning (disable: 4996)
//#pragma warning (disable: LNK4229)

#define ESCape 0x1b

void VVOD ( void );
void F_df(int q);
void DROITE(void);
void D_N_PER( void );
void ZERO( double **j,double **J);
void CONVERG( void);
void MATR_(  int mz,int mz1);
//void MATR_S( void);
void C_sqrt(double x,double y, double *u,double *v);
void C_exp(double x,double y, double *u,double *v);
void C_sh(double x,double y, double *u,double *v);
void C_ch(double x,double y, double *u,double *v);
void C_sh_ch(double x,double y, double *us,double *vs, double *uc,double *vc);
void Cn_sh_ch(double rg_,double ig_,double z,double h, double *us,double *vs, double *uc,double *vc);
void EXTERN(double ****E, double ****H);
void sq_compl(double x,double y,double *u,double *v);
void SH_CTH(double x,double y,double *us,double *vs,double *uc,double *vc);
void SH_CTH_2(double ro,double x,double y,double *us,double *vs,double *us2,double *vs2,double *uc,double *vc);
void sin_cos_C(double x,double y,double *us,double *vs,double *uc,double *vc);
void sin_C(double x,double y,double *us,double *vs);
void sinh_C(double x,double y,double *us,double *vs);
void cosh_C(double x,double y,double *us,double *vs);
void CONST(int N_d,int i,double Ro, double *te,double *qu,double *fi1,double *fi3,double *del,double *S1,double *S3);
double Max(double x, double y);
double Ag(double x, int ri);
double Au(double x, int ri);
double Cu(double x, int ri);
void DIAGON( void);
void DIAGON_2( void);
void NOMER( void);
void E_extern ( double *kx_,double *ky_,double *kz_  );
void Ancillary( void);
void INT(int N_d,int h,int m, double kap_r, double kap_i, double **S, double **C, double *s, double *c);
void BASIS( double alf, double bet1);//new
void ZERO_B(void);//new
void FM1( void);
void AM1( void );
void FM2( void);
void AM2( void );
void FM3( void);
void AM3( void );
void FM4( void);
void AM4( void );
void FM5( void);
void AM5( void );
void FM6( void);
void AM6( void );
void FM0( void);
void AM0( void );
void DIAGON_N( void);
void RECURS(double Ro,int N_d,int i,double **rte,double **ite,double **rqu,double **iqu, double *rG1,double *iG1, double *rG3,double *iG3);
void RECURS_ZERO(double Ro,int N_d, int i, double **rte, double **ite, double **rqu, double **iqu, double *rG1, double *iG1, double *rG3, double *iG3);
void COMPOSITE_me(double R_eps, double I_eps, double C, int n_mat, double *R_eps_c, double *I_eps_c);
void COMPOSITE(double R_eps, double I_eps, double C, int n_mat, double *R_eps_c, double *I_eps_c);

void COMPOSITE_3(double R_eps_m_1, double R_eps_m_2, double I_eps_m_1, double I_eps_m_2, double C1, double C2, int n_mat_1, int n_mat_2, double* R_eps_c, double* I_eps_c);
/*
double Ag(double x, int ri);
double Au(double x, int ri);
double Cu(double x, int ri);
void ZnO_(double lam, double *Re_eps, double *Im_eps);
void Ag_(double lam, double *Re_eps, double *Im_eps);
void Au_(double lam, double *Re_eps, double *Im_eps);
void Si_(double lam, double *Re_eps, double *Im_eps);
void SiO2_(double lam, double *Re_eps, double *Im_eps);
void HfO_(double lam, double *Re_eps, double *Im_eps);
void AZO_(double lam, double *Re_eps, double *Im_eps);
double Ag(double x, int ri);
double Au(double x, int ri);
double Cu(double x, int ri);
void ZnO_(double lam, double *Re_eps, double *Im_eps);
void Ag_(double lam, double *Re_eps, double *Im_eps);
void Au_(double lam, double *Re_eps, double *Im_eps);
void Si_(double lam, double *Re_eps, double *Im_eps);
void SiO2_(double lam, double *Re_eps, double *Im_eps);
void HfO_(double lam, double *Re_eps, double *Im_eps);
void AZO_(double lam, double *Re_eps, double *Im_eps);

void VO2_cold(double x, double *eps1, double *eps2);
void VO2_hot(double x, double *eps1, double *eps2);
*/
__inline void  mult_( double x1,double y1,double x2,double y2,double *u,double *v )
{
	*u = x1 * x2 - y1 * y2;
	*v = x1 * y2 + y1 * x2;
}

__inline void divide_(double x,double y,double x1,double y1,double *u,double *v)
{
	double z;
	z = 1.0 / ( x1 * x1 + y1 * y1 );
	*u = ( x * x1 + y * y1 ) * z;
	*v = ( -x * y1 + y * x1 ) * z;
}
//void PRO_SIM( void);

int *m_t,*m_b,ms,M,AM[4],graf1,graf2,d_n,Nv,MS[2],MS_,MS1,NS,N_corn,DN,Pol,M_z,M_fi,M_r,AM_z[4],AM_fi[4],AM_r[4],
N_g, AN_g[4], ****No, ****No2, mfi1[2], mfi2[2], M_[4], e, SIM_r, *N_r, Mz, *is, *comp_,is_in, *N_l, Nsl, N_hol, bur_v, PA,
	bu_p[201],Ng_r,BU_SIM,
	M_r_max,N_La,M_r_,**N_lin,*MR,M_r_1, Ekv,N_az,N_rad,
	N_ell, MS_N, is_l[21][6], comp_l[21][6], bum, bum1;

double k,k2,ky,kz,kx,AF[4],f,Pi_2,Pi,TET,ATET[4],AFI[4],teta,fi,FI,FI_r,E0[3][4],H0[3][4],****E,****H,****Ep,****Hp,
	   Z0,**ASI_,**ASR_,Ce[3],*a_,*b_,t_,DE1[3],DEN1[3],*n_,*k_,*H_,n_in,k_in,EPS_r_in,EPS_i_in,t_r,tet0,
	   *h_, *h_b, *EPS_r, *EPS_i, *EPS_r_, *EPS_i_, dx, dy, k0, cotet, sitet, cofi, sifi, TET_r, rX[5], iX[5], FI_r_,
		***R_c2,***R_s2,**g_r,**g_i,*k2_r,*k2_i,*k_r,*k_i,SQ_05,Z0_,lam,a_in,b_in,H_s,
		kx,ky,kz,pi2,***j,***J,h,Ex,Ey,Ez,*A_t,*B_t,*A_b,*B_b,*I_s,*I_c,*I_s_i,*I_c_i,
		Ex_0, Ey_0, Ez_0, Hx_0, Hy_0, Hz_0, *CONCEN, CONCEN_l[21][6],
		*Zm,*hm,*hm1,*dm,*rg_,*ig_,**Gm1,**Sh,**Ch,***U,***V,***Wu,***Wv,**Vd, ****D_t,****D_b,
		**A_t_l, **B_t_l, **A_b_l, **B_b_l, **EPS_l_r, **EPS_l_i, **EPS_l_r_, **EPS_l_i_, **n_l, **k_l, ft_l,
		
		 **er,**ei,G_r[10],G_i[10],Rad_ro,
		
		**r_l,*fi_0,*c_fi0,*s_fi0,***Y1,***A_k,**al_,**bl_,*XC,*YC,*Co_c,*Si_c, dze[4][111], dze_im[4][111],tg_pci;

scomplex **Alap,**Blap;


FILE *in,*out,*out1,*out2,*in1;
char inf[8], outf[8], out1f[8], out2f[8], in1f[8];
clock_t str,start;
void tim ( clock_t str );

