
#include "VIBRr.h"

void AM0( void)//new
{
	E = ( double****)amat_4d( 2,N_hol+1,1,2,1,3, 1,4, sizeof( double) );
	H = ( double****)amat_4d( 2,N_hol+1,1,2,1,3, 1,4, sizeof( double) );
	Ep = ( double****)amat_4d( 2,N_hol+1,1,2,1,3, 1,4, sizeof( double) );
	Hp = ( double****)amat_4d( 2,N_hol+1,1,2,1,3, 1,4, sizeof( double) );
}

void FM0( void)//new
{
	fmat_4d( (void****)E, 2,N_hol+1,1,2,1,3, 1, sizeof( double) );
	fmat_4d( (void****)H, 2,N_hol+1,1,2,1,3, 1, sizeof( double) );
	fmat_4d( (void****)Ep, 2,N_hol+1,1,2,1,3, 1, sizeof( double) );
	fmat_4d( (void****)Hp, 2,N_hol+1,1,2,1,3, 1, sizeof( double) );
}

void AM1( void)//new
{
	n_ = ( double*)amat_1d( 0,Nsl + 1, sizeof( double) );
	k_ = ( double*)amat_1d( 0,Nsl + 1, sizeof( double) );
	H_ = ( double*)amat_1d( 0,Nsl + 1, sizeof( double) );
	EPS_r = ( double*)amat_1d( 0,Nsl + 1, sizeof( double) );
	EPS_i = ( double*)amat_1d( 0,Nsl + 1, sizeof( double) );
	EPS_r_ = (double*)amat_1d(0, Nsl + 1, sizeof(double));
	EPS_i_ = (double*)amat_1d(0, Nsl + 1, sizeof(double));
	k2_r = ( double*)amat_1d( 0,Nsl + 1, sizeof( double) );
	k2_i = ( double*)amat_1d( 0,Nsl + 1, sizeof( double) );
	k_r = ( double*)amat_1d( 0,Nsl + 1, sizeof( double) );
	k_i = ( double*)amat_1d( 0,Nsl + 1, sizeof( double) );
	rg_ = ( double*)amat_1d( 0,Nsl + 1, sizeof( double) );
	ig_ = ( double*)amat_1d( 0,Nsl + 1, sizeof( double) );
	h_ = ( double*)amat_1d( 0,Nsl + 1, sizeof( double) );
	h_b = ( double*)amat_1d( 0,Nsl + 1, sizeof( double) );

	is = ( int*)amat_1d( 0,Nsl + 1, sizeof( int) );
	comp_ = (int*)amat_1d(0, Nsl + 1, sizeof(int));
	
	CONCEN = (double*)amat_1d(0, Nsl + 1, sizeof(double));

	Gm1 = ( double**)amat_2d(1,2,0,Nsl + 1, sizeof( double) );
}

void FM1( void)//new
{
	fmat_1d( (void*)n_,  0, sizeof( double) );
	fmat_1d( (void*)k_,  0, sizeof( double) );
	fmat_1d( (void*)H_,  0, sizeof( double) );
	fmat_1d( (void*)EPS_r,  0, sizeof( double) );
	fmat_1d( (void*)EPS_i,  0, sizeof( double) );
	fmat_1d((void*)EPS_r_, 0, sizeof(double));
	fmat_1d((void*)EPS_i_, 0, sizeof(double));
	fmat_1d( (void*)k_r,  0, sizeof( double) );
	fmat_1d( (void*)k_i,  0, sizeof( double) );
	fmat_1d( (void*)k2_r,  0, sizeof( double) );
	fmat_1d( (void*)k2_i,  0, sizeof( double) );
	fmat_1d( (void*)rg_,  0, sizeof( double) );
	fmat_1d( (void*)ig_,  0, sizeof( double) );
	fmat_1d( (void*)h_,  0, sizeof( double) );
	fmat_1d( (void*)h_b,  0, sizeof( double) );
	
	fmat_1d( (void*)is,  0, sizeof( int) );
	fmat_1d((void*)comp_, 0, sizeof(int));
	fmat_1d((void*)CONCEN, 0, sizeof(double));

	fmat_2d( (void**)Gm1, 1,2, 0, sizeof( double) );
}

void AM2( void)//new
{
	A_t = ( double*)amat_1d( 1,N_hol, sizeof( double) );
	B_t = ( double*)amat_1d( 1,N_hol, sizeof( double) );
	A_b = ( double*)amat_1d( 1,N_hol, sizeof( double) );
	B_b = ( double*)amat_1d( 1,N_hol, sizeof( double) );
	fi_0 = ( double*)amat_1d( 1,N_hol, sizeof( double) );
	c_fi0 = ( double*)amat_1d( 1,N_hol, sizeof( double) );
	s_fi0 = ( double*)amat_1d( 1,N_hol, sizeof( double) );
	
	m_t = ( int*)amat_1d( 1,N_hol+1, sizeof( int) );
	m_b = ( int*)amat_1d( 1,N_hol+1, sizeof( int) );
	N_l = ( int*)amat_1d( 1,N_hol, sizeof( int) );
	MR = ( int*)amat_1d( 1,N_hol, sizeof( int) );
}

void FM2( void)//new
{
	fmat_1d( (void*)A_t,  1, sizeof( double) );
	fmat_1d( (void*)B_t,  1, sizeof( double) );
	fmat_1d( (void*)A_b,  1, sizeof( double) );
	fmat_1d( (void*)B_b,  1, sizeof( double) );
	fmat_1d( (void*)fi_0,  1, sizeof( double) );
	fmat_1d( (void*)c_fi0,  1, sizeof( double) );
	fmat_1d( (void*)s_fi0,  1, sizeof( double) );
	
	fmat_1d( (void*)m_t,  1, sizeof( int) );
	fmat_1d( (void*)m_b,  1, sizeof( int) );
	fmat_1d( (void*)N_l,  1, sizeof( int) );
	fmat_1d( (void*)MR,  1, sizeof( int) );
}

void AM3( void)//new
{
	A_t_l = ( double**)amat_2d(1,N_hol,1,N_La, sizeof( double) );
	B_t_l = ( double**)amat_2d(1,N_hol,1,N_La, sizeof( double) );
	A_b_l = ( double**)amat_2d(1,N_hol,1,N_La, sizeof( double) );
	B_b_l = ( double**)amat_2d(1,N_hol,1,N_La, sizeof( double) );
	EPS_l_r = ( double**)amat_2d(1,N_hol,1,N_La, sizeof( double) );
	EPS_l_i = ( double**)amat_2d(1,N_hol,1,N_La, sizeof( double) );
	EPS_l_r_ = (double**)amat_2d(1, N_hol, 1, N_La, sizeof(double));
	EPS_l_i_ = (double**)amat_2d(1, N_hol, 1, N_La, sizeof(double));
	n_l = ( double**)amat_2d(1,N_hol,1,N_La, sizeof( double) );
	k_l = ( double**)amat_2d(1,N_hol,1,N_La, sizeof( double) );
	XC = ( double*)amat_1d(1,N_ell, sizeof( double) );
	YC = ( double*)amat_1d(1,N_ell, sizeof( double) );
	Co_c = ( double*)amat_1d(1,N_ell, sizeof( double) );
	Si_c = ( double*)amat_1d(1,N_ell, sizeof( double) );
	er = ( double**)amat_2d(1,N_ell,1,N_ell, sizeof( double) );
	ei = ( double**)amat_2d(1,N_ell,1,N_ell, sizeof( double) );
}


void FM3( void)//new
{
	fmat_2d( (void**)A_t_l, 1,N_hol, 1, sizeof( double) );
	fmat_2d( (void**)B_t_l, 1,N_hol, 1, sizeof( double) );
	fmat_2d( (void**)A_b_l, 1,N_hol, 1, sizeof( double) );
	fmat_2d( (void**)B_b_l, 1,N_hol, 1, sizeof( double) );
	fmat_2d( (void**)EPS_l_r, 1,N_hol, 1, sizeof( double) );
	fmat_2d( (void**)EPS_l_i, 1,N_hol, 1, sizeof( double) );
	fmat_2d((void**)EPS_l_r_, 1, N_hol, 1, sizeof(double));
	fmat_2d((void**)EPS_l_i_, 1, N_hol, 1, sizeof(double));
	fmat_2d( (void**)n_l, 1,N_hol, 1, sizeof( double) );
	fmat_2d( (void**)k_l, 1,N_hol, 1, sizeof( double) );
	fmat_1d( (void*)XC,  1, sizeof( double) );
	fmat_1d( (void*)YC,  1, sizeof( double) );
	fmat_1d( (void*)Co_c,  1, sizeof( double) );
	fmat_1d( (void*)Si_c,  1, sizeof( double) );
	fmat_2d( (void**)er, 1,N_ell, 1, sizeof( double) );
	fmat_2d( (void**)ei, 1,N_ell, 1, sizeof( double) );
}

void AM4( void)//new
{
	N_lin = ( int**)amat_2d( 1,N_hol,1,M_r_max, sizeof( int) );
	j = ( double***)amat_3d( 0,AM_fi[3]-1, 1,M_r_max, 1,N_hol, sizeof( double) );
	J = ( double***)amat_3d(  0,AM_fi[3]-1, 1,M_r_max, 1,N_hol,sizeof( double) );
	if (N_La > 1)
	{
		Y1 = ( double***)amat_3d(  0,AM_fi[3]-1, 1,M_r_max, 1,N_hol,sizeof( double) );
		A_k = ( double***)amat_3d(  0,AM_fi[3]-1, 1,M_r_max, 1,N_hol,sizeof( double) );
	}
}

void FM4( void)//new
{
	fmat_2d( (void**)N_lin,1,N_hol,  1, sizeof( int) );
	fmat_3d( (void***)j, 0,AM_fi[3]-1, 1,M_r_max, 1, sizeof( double) );
	fmat_3d( (void***)J, 0,AM_fi[3]-1, 1,M_r_max, 1, sizeof( double) );
	if (N_La > 1)
	{
		fmat_3d( (void***)Y1, 0,AM_fi[3]-1, 1,M_r_max, 1, sizeof( double) );
		fmat_3d( (void***)A_k, 0,AM_fi[3]-1, 1,M_r_max, 1, sizeof( double) );
	}
}

void AM5( void)//new
{
	ASR_ = ( double**)amat_2d(0,1,1,MS_N, sizeof( double) );// third party code (source belongs to A.L. colleague)
	ASI_ = ( double**)amat_2d(0,1,1,MS_N, sizeof( double) );
	Zm = ( double*)amat_1d( 0,ms+1, sizeof( double) );
	hm = ( double*)amat_1d( 1,ms+1, sizeof( double) );
	hm1 = ( double*)amat_1d( 1,ms+1, sizeof( double) );
	dm = ( double*)amat_1d( 1,ms, sizeof( double) );
	a_ = ( double*)amat_1d( 1,ms, sizeof( double) );
	b_ = ( double*)amat_1d( 1,ms, sizeof( double) );

	al_ = ( double**)amat_2d( 1,ms, 1,N_La, sizeof( double) );
	bl_ = ( double**)amat_2d( 1,ms, 1,N_La, sizeof( double) );

	r_l = ( double**)amat_2d( 1,N_hol,0,N_La, sizeof( double) );//new

	N_r = ( int*)amat_1d( 0,ms+1, sizeof( int) );
	// memory allocation for 4D-array
	No = ( int****)amat_4d( 1,3,1,ms, 1,M_r, 1,M_fi,  sizeof( int) );
	No2 = ( int****)amat_4d( 1,3, 1,ms,1,M_r, 1,M_fi,  sizeof( int) );
}

void FM5( void)//new
{
	DeAllocC(Alap);
	DeAllocC(Blap);
			
	fmat_1d( (void*)N_r,  0, sizeof( int) );
	fmat_1d( (void*)a_,  1, sizeof( double) );
	fmat_1d( (void*)b_,  1, sizeof( double) );

	fmat_2d( (void**)al_, 1,ms, 1, sizeof( double) );
	fmat_2d( (void**)bl_, 1,ms, 1, sizeof( double) );

	fmat_2d( (void**)r_l,1,N_hol,0, sizeof( double) );
	fmat_2d( (void**)ASI_, 0,1, 1, sizeof( double) );//
	fmat_2d( (void**)ASR_, 0,1, 1, sizeof( double) );//
	fmat_1d( (void*)Zm,  0, sizeof( double) );
	fmat_1d( (void*)hm,  1, sizeof( double) );
	fmat_1d( (void*)hm1,  1, sizeof( double) );
	fmat_1d( (void*)dm,  1, sizeof( double) );
	fmat_4d( (void****)No,1,3, 1,ms,1,M_r, 1, sizeof( int) );
	fmat_4d( (void****)No2,1,3,1,ms,1,M_r, 1, sizeof( int) );
}


void AM6( void)//new
{
	 R_c2 = ( double***)amat_3d( 1,M_z,1,M_fi, 1,M_r, sizeof( double) );
	R_s2 = ( double***)amat_3d( 1,M_z,1,M_fi, 1,M_r, sizeof( double) );
	
	g_r = ( double**)amat_2d( 1,9,1,M_z,  sizeof( double) );
	g_i = ( double**)amat_2d( 1,9,1,M_z,  sizeof( double) );

	I_s = ( double*)amat_1d( 1,M_z, sizeof( double) );
	I_c = ( double*)amat_1d( 1,M_z, sizeof( double) );
	I_s_i = ( double*)amat_1d( 1,M_z, sizeof( double) );
	I_c_i = ( double*)amat_1d( 1,M_z, sizeof( double) );
	
	Sh = ( double**)amat_2d( 1,2, 0,M_z+1, sizeof( double) );
	Ch = ( double**)amat_2d( 1,2, 0,M_z+1, sizeof( double) );
	
	
	U = ( double***)amat_3d( 1,2,0,M_z+1, 0,M_z+1, sizeof( double) );
	V = ( double***)amat_3d( 1,2,0,M_z+1, 0,M_z+1, sizeof( double) );
	Wu = ( double***)amat_3d( 1,2,0,M_z+1, 0,M_z+1, sizeof( double) );
	Wv = ( double***)amat_3d( 1,2,0,M_z+1, 0,M_z+1, sizeof( double) );
	Vd = ( double**)amat_2d( 1,2, 0,ms+1, sizeof( double) );

	D_t = ( double****)amat_4d( 1,2,1,3, 1,N_hol+1, 1,N_hol+1, sizeof( double) );
	D_b = ( double****)amat_4d( 1,2,1,3, 1,N_hol+1, 1,N_hol+1, sizeof( double) );
}

void FM6( void)//new
{
	fmat_2d( (void**)g_r,1,9,  1, sizeof( double) );
	fmat_2d( (void**)g_i,1,9,  1, sizeof( double) );

	fmat_3d( (void***)R_s2,1,M_z, 1,M_fi, 1, sizeof( double) );
	fmat_3d( (void***)R_c2,1,M_z, 1,M_fi, 1, sizeof( double) );
	
	//fmat_1d( (void*)J_n,  0, sizeof( double) );//new
	fmat_1d( (void*)I_s,  1, sizeof( double) );
	fmat_1d( (void*)I_c,  1, sizeof( double) );
	fmat_1d( (void*)I_s_i,  1, sizeof( double) );
	fmat_1d( (void*)I_c_i,  1, sizeof( double) );
	fmat_2d( (void**)Sh, 1,2, 0, sizeof( double) );
	fmat_2d( (void**)Ch, 1,2, 0, sizeof( double) );

	fmat_3d( (void***)U,1,2, 0,M_z+1, 0, sizeof( double) );
	fmat_3d( (void***)V,1,2, 0,M_z+1, 0, sizeof( double) );
	fmat_3d( (void***)Wu,1,2, 0,M_z+1, 0, sizeof( double) );
	fmat_3d( (void***)Wv,1,2, 0,M_z+1, 0, sizeof( double) );
	fmat_2d( (void**)Vd,1,2, 0, sizeof( double) );

	fmat_4d( (void****)D_t, 1,2,1,3, 1,N_hol+1,1, sizeof( double) );
	fmat_4d( (void****)D_b, 1,2,1,3, 1,N_hol+1,1, sizeof( double) );
}

