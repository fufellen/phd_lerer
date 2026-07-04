
#include "VIBRr.h"

//double Ir(int p,int m1,int m2,double r);
void NORMA(  int n_hol,int mf, int mr , int mr1, double R1,double R2,double *r_N  );
int M;


void DIAGON( void)//Диагональные элементы.   Аналитическое  интегрирование.
{
	int N_d,n,i,mz1,mz,mf,mr,mr1,m,m1,c,mb,mt,P,nl,
		n_hol;
	double x,y,z,z1,z2,u_,v_,R1;
	double **u,**v;
	u = ( double**)amat_2d(1,N_hol +1,1,N_La, sizeof( double) );
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
	{
		for ( mz=1; mz <= M_z; mz++ )//new
		{
			N_d = N_r[mz];//Номер слоя, в котором расположена точка Zm[m1]

			n_hol = N_d - 1;

			mb = m_b[N_d];
			mt = m_t[N_d];
				
			
			for ( mr=1; mr <= MR[n_hol]; mr++ )//new
			//for ( mr1=1; mr1 <= M_r; mr1++ )
			{
				mr1 = mr;
				nl = N_lin[n_hol][mr];	//new
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
					NORMA(N_d-1,mf, mr , mr1, R1,r_l[N_d-1][nl],&z);
					
					u_ = u[N_d][nl] * z;
					v_ = v[N_d][nl] * z;
	
					if (mf == 1)
					{
						u_ *= 2.0;
						v_ *= 2.0;
					}

					for ( P=0; P <= 1; P++ )
					for ( i=1; i <= 3; i++ )
					{
						c = (i == 1) ? 1 + P : 2 - P;
						if (mf >= c)
						{
							m = (P) ? No2[i][mz][mr][mf] : No[i][mz][mr][mf];
							m1 = (P) ? No2[i][mz1][mr1][mf] : No[i][mz1][mr1][mf];
							
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
}//void DIAGON( void)

void NORMA(  int n_hol, int mf, int mr , int mr1, double R1,double R2,double *r_N  )//норма  функций Бесселя
{
	double b,b1,U,v;
	
	b = J[mf-1][mr][n_hol] * R2;
	b1 = (mf == 1) ? 0.0 : (double)(mf-1) / j[mf-1][mr][n_hol];
	U = b1 / R2;
	U = 0.5 * b * b * (1.0 - U * U); 
	
	if (mr > M)
	{
		b = Y1[mf-1][mr][n_hol] * R1;
		v = b1 / R1;
		U -= 0.5 * b * b * (1.0 - v * v); 
	}
	*r_N = U;
	
}//Norma
	
