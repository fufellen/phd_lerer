
#include "VIBRr.h"

void NORMA_N(  int P,int i, int n_hol,int nl,int mf, int mf1, int mr , int mr1, int mz, double R1,double R2,double *u_,double *v_  );
int M;
double **u,**v;

void DIAGON_N( void)//Диагональные элементы. Численное интегрирование.
{
	int N_d,n,i,mz1,mz,mf,mf1,mr,mr1,m,m1,c,mb,mt,P,
		n_hol,nl,nl1;
	double x,y,z1,z2,u_,v_,R1,R2;
		
	u = ( double**)amat_2d(1,N_hol+1,1,N_La, sizeof( double) );
	v = ( double**)amat_2d(1,N_hol+1,1,N_La, sizeof( double) );

	//_______________________________________________________________________
	z2 = -dx * dy / (8.0 * Pi);

	for ( n=1; n <= N_hol; n++ )
		for ( m=1;m <= N_l[n]; m++ )
		{
			N_d = n + 1;
			divide_(z2,0.0, EPS_l_r[n][m] - EPS_r[N_d],EPS_l_i[n][m] - EPS_i[N_d], &x,&y);
			mult_(x,y, EPS_r[N_d],EPS_i[N_d], &u[N_d][m],&v[N_d][m]);
		}
	
	for ( mf=1; mf <= M_fi; mf++ )
	for ( mf1=1; mf1 <= M_fi; mf1++ )
	{
		//mf1 = mf;
		for ( mz=1; mz <= M_z; mz++ )//new
		{
			N_d = N_r[mz];//Номер слоя, в котором расположена точка Zm[m1]

			n_hol = N_d - 1;

			mb = m_b[N_d];
			mt = m_t[N_d];
				
			
			for ( mr=1; mr <= MR[n_hol]; mr++ )//new
			for ( mr1=1; mr1 <= MR[n_hol]; mr1++ )
			{
				//mr1 = mr;
				nl = N_lin[n_hol][mr];	
				nl1 = N_lin[n_hol][mr1];
				if (nl != nl1)
					continue;

				M = (nl == 1) ? M_r_1 : M_r_;//new
													
				for ( mz1=1; mz1 <= M_z; mz1++ )
				{
					if (N_d != N_r[mz1])
						continue;

					if (mz1 == mz)
						z1 = (mz == mt) ? hm[mz] / 3.0 
											: (mz == mb) ? hm[mz+1] / 3.0
														: (hm[mz] + hm[mz+1]) / 3.0;
						
					else
						z1 =  (mz1 == mz + 1) ? hm[mz+1] / 6.0
								    		    : (mz1 == mz - 1) ? hm[mz] / 6.0 : 0.0;
					if (z1 < 0.00001)
						continue;
						
					z1 /= a_[mz1] * b_[mz1];
					
					R1 = (nl > 1) ? r_l[N_d-1][nl-1] : 0.0;
					R2 = r_l[N_d-1][nl];
										
					for ( P=0; P <= 1; P++ )
					for ( i=1; i <= 3; i++ )
					{
						c = (i == 1) ? 1 + P : 2 - P;

						if (mf >= c)
						{
							NORMA_N(P,i,N_d-1,nl,mf,mf1, mr , mr1, mz1,R1,R2,&u_,&v_);
							
							m = (P) ? No2[i][mz][mr][mf] : No[i][mz][mr][mf];
							m1 = (P) ? No2[i][mz1][mr1][mf1] : No[i][mz1][mr1][mf1];
						
	
							Alap[m][m1].real += u_ * z1;
							Alap[m][m1].imag += v_ * z1;
						}
					}//i,P
						
				}//mz1
					
			}//mz
		}//mr
	}//mf
	
	fmat_2d( (void**)u, 1,N_hol+1, 1, sizeof( double) );
	fmat_2d( (void**)v, 1,N_hol+1, 1, sizeof( double) );

}//DIAGON_N

void NORMA_N(  int P,int I, int n_hol,int nl,int mf,int mf1, int mr , int mr1, int mz, double R1,double R2,double *u_,double *v_ )
{
	int i,l,N_d,nl_i,bu_s,N_,N_n;
	double x,y,ro,fi,s,z,b,b1,h,h_fi,r,U,V,mf_,mf1_,a,be, U_,V_,Uu,Vv, Uu1,Vv1,si[1020],co[1020],si_m[1020],si1_m[1020],
		z0,z01,A,A1;
	
	N_= N_az;
	N_n = N_rad;

	N_d = n_hol + 1;
	a = a_[mz] / al_[mz][nl];
	be = b_[mz] / bl_[mz][nl];

	z0 = j[mf-1][mr][n_hol];
	z01 = j[mf1-1][mr1][n_hol];
	A = A_k[mf-1][mr][n_hol];
	A1 = A_k[mf1-1][mr1][n_hol];
	
	Uu = u[N_d][nl];
	Vv = v[N_d][nl];
	Uu1 = u[N_d][nl+1];
	Vv1 = v[N_d][nl+1];

	mf_ = (double)(mf - 1);
	mf1_ = (double)(mf1 - 1);

	bu_s = ( (I == 1 && P) || ((I != 1 && !P)) ) ? 1 : 0;
		
	h_fi = 2.0 * Pi / N_;

	for ( l=1; l <= N_; l++ )
	{
		fi = (-0.5 + l) * h_fi;
			
		co[l] = cos(fi);
		si[l] = sin(fi);
			
		if (bu_s)
		{
			si_m[l] = sin(fi * mf_);
			si1_m[l] = sin(fi * mf1_);
		}
		else
		{
			si_m[l] = cos(fi * mf_);
			si1_m[l] = cos(fi * mf1_);
		}
	}//l

	h = (R2 - R1) / N_n;
	for ( V=U=0, i=1; i <= N_n; i++ )
	{
		r = h * (-0.5 + i) + R1;
	
		b = (z0 < 1.0e-17) ? 1.0 : (mr <= M) ? jn(mf-1,z0 * r)  : jn(mf-1,z0 * r) +  A * yn(mf-1,z0 * r);
		b1 = (z01 < 1.0e-17) ? 1.0 : (mr1 <= M) ? jn(mf1-1,z01 * r)  : jn(mf1-1,z01 * r) +  A1 * yn(mf1-1,z01 * r);

		z = r * b * b1;
				
		for ( l=1; l <= N_; l++ )
		{
			fi = (-0.5 + l) * h_fi;
			
			x = a * r * co[l];
			y = be * r * si[l];
			ro = x * x + y * y;

			if (ro <= 1.0)
			{
				nl_i = nl;
				U_ = Uu;
				V_ = Vv;
			}
			else
			{
				nl_i = nl + 1;
				U_ = Uu1;
				V_ = Vv1;
			
			}
								
			s = z * si_m[l] * si1_m[l];

			U += U_ * s;
			V += V_ * s;
		}//l
	}//i
		
	z = h * h_fi / Pi;
	 
	*u_ = U * z;
	*v_ = V * z;
	
}//Norma
