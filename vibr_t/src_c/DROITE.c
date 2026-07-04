
#include "VIBRr.h"

void DROITE(void)//внешнее поле exp(i(kx*x+ky*y+kz*z)).Амплитуда падающей волны 1 В / м external field is exp (i (kx * x + ky * y + kz * z)). The amplitude of the incident wave is 1 V / m

{ // calc. Integral{Ee*V(mfi,mr,mz)}dV for external field Ee multiplied by Base function over elliptic cylinder volume.
	// the integral is taken analytically. External field is a plane wave "zero"-harmonic.
	// result is in Blap i.e. right side of the equation system;
	// Blap is 2D array {illumIndex, runIndexOver [(X,Y,Z)l,m,n] }
	// Blap can be 1D because we solved the system for a single illumination direction
	double n1,u,v,u1,v1,x,y,z,yy,sz[3],cz[3],sh[3],ch[3],A[3],B[3],**S,**CC,xc,yc,kap_r[221],kap_i[221];//new
	int N_d,c,i,mf,n,m,mr,mz,mz_,P,
		n_hol,nl,m_;

	S = ( double**)amat_2d( 1,2, 0,ms+1, sizeof( double) );
	CC = ( double**)amat_2d( 1,2, 0,ms+1, sizeof( double) );

		
	for ( mz=1; mz <= ms; mz++ )
	{
		N_d = N_r[mz];
		xc =  k2_r[N_d] - ky * ky - kx * kx;
		yc =  k2_i[N_d];
		C_sqrt(xc,yc,&kap_r[N_d],&kap_i[N_d]);

		x = (Zm[mz] - h_b[N_d]) * kap_r[N_d];
		y = (Zm[mz] - h_b[N_d]) * kap_i[N_d];
		sin_cos_C(x,y, &u,&v, &u1,&v1);
		S[1][mz] = u;
		S[2][mz] = v;
		CC[1][mz] = u1;
		CC[2][mz] = v1;
	}

	BASIS(kx,ky);//вычисление Rnp calculate

	for ( nl=1; nl <= N_ell; nl++ )
	{
		u = kx * XC[nl] + ky * YC[nl];
		Co_c[nl] = cos(u);
		Si_c[nl] = sin(u);
	}

	for ( Pol=0; Pol <= 1; Pol++ )
	for ( P=0; P <= 1; P++ )
	{
		for ( i=1; i <= 2 + Pol; i++ )
		{
			c = (i == 1) ? 1 + P : 2 - P;
			
			for ( mf=c; mf <= M_fi; mf++ )
			{
				n = mf - 1;
				n1 = (double)n;		
			
				for ( mz=1; mz <= ms; mz++ )//new
				{
					n_hol = N_r[mz]-1;
					
					for ( mr=1; mr <= MR[n_hol]; mr++ )//new
					{
						N_d = N_r[mz];

						if (Pol)
						{
							A[1] = -Ep[N_d][1][i][3];
							A[2] = -Ep[N_d][2][i][3];
							B[1] = -Ep[N_d][1][i][4];
							B[2] = -Ep[N_d][2][i][4];
						}
						else
						{
							A[1] = -E[N_d][1][i][3];
							A[2] = -E[N_d][2][i][3];
							B[1] = -E[N_d][1][i][4];
							B[2] = -E[N_d][2][i][4];
							
						}
					
										
						INT(N_d,0,mz, kap_r[N_d], kap_i[N_d], S, CC, sz, cz);

						mz_ = m_t[N_d] - (mz - m_b[N_d]);
				//		if (N_d != N_r[mz_])
				//			printf ("\n OSHIBKA 4\n");
					
						INT(N_d,0,mz_, kap_r[N_d], kap_i[N_d], S, CC, sh, ch);
					
						//new
						if (P)
							yy = (i == 1)? R_s2[mz][mf][mr] : R_c2[mz][mf][mr];
						else
							//yy = (i == 1)? R_s[mz][mf][mr] : R_c[mz][mf][mr];
							yy = (i == 1)? R_c2[mz][mf][mr] : R_s2[mz][mf][mr];//new
									
						//new end
						
						if (Pol)
						{
							if (i == 3)
							{
								mult_(A[1],A[2], sz[1],sz[2], &u,&v);
								mult_(B[1],B[2], sh[1],sh[2], &u1,&v1);
								u += u1;
								v += v1;
							}
							else
							{
								mult_(A[1],A[2], cz[1],cz[2], &u,&v);
								mult_(B[1],B[2], ch[1],ch[2], &u1,&v1);
								u -= u1;
								v -= v1;
							}
							if (i == 3)
							{
								z = u;
								u = -v;
								v = z;
							}
						}//Pol
						else
						{
								mult_(A[1],A[2], sz[1],sz[2], &u,&v);
								mult_(B[1],B[2], sh[1],sh[2], &u1,&v1);
								u += u1;
								v += v1;
						}
						
						m = (P) ? No2[i][mz][mr][mf] : No[i][mz][mr][mf];
										
						Blap[Pol][m].real = yy * u;
						Blap[Pol][m].imag = yy * v;
					
					}//mr
				}//mz
			}//mf
		}//i
	}//P

	for ( Pol=0; Pol <= 1; Pol++ )
		for ( m=0; m < MS_; m++ )
		{
			u = Blap[Pol][m].real;
			v = Blap[Pol][m].imag;

			for ( i=1; i <= N_ell; i++ )//неоднородная часть The inhomogeneous part
			{
				m_ = m + MS_* (i - 1);
				
				mult_(u,v, Co_c[i], Si_c[i], &u1,&v1);
				Blap[Pol][m_].real = u1;
				Blap[Pol][m_].imag = v1;

			}//i
		}//m
			
	fmat_2d( (void**)S, 1,2, 0, sizeof( double) );
	fmat_2d( (void**)CC, 1,2, 0, sizeof( double) );
}/* DROITE */

