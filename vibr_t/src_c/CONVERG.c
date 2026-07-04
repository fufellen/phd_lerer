
#include "VIBRr.h"

double delta_K(int m);
double Ii(int P,int i,int i1,int m,int mz);
double R_m;

void CONVERG( void)
{ // asymptotic calc. of integrals outside diffr.scheme
	int N_d,mb,mt,c,c1,i,i1,mz,mz1,mf,mr,mr1,n,m,m1,P,
		n_hol,nl,nl1,M;
	double AL,u,x,z,u2,
		B,Rv,Rv1,Ry,Ry1;
			
	
	////R_m = Max(Pi_2 * (double)N_g / dx,Pi_2 * (double)N_g / dy);
	R_m = (Pi_2 * (double)N_g / dx + kx < Pi_2 * (double)N_g / dy + ky) ? Pi_2 * (double)N_g / dx  + kx 
																										        : Pi_2 * (double)N_g / dy + ky;

	
	//___________________________________________________________________

	for ( P=0; P <= 1; P++ )
	for ( i=1; i <= 3; i++ )//улучшение сходимости improvement of convergence
	{
		c = (i == 1) ? 1 + P : 1 + ( 1 - P);
	
		for ( i1=1; i1 <= 3; i1++ )//улучшение сходимости improvement of convergence
		{
			c1 = (i1 == 1) ? 1 + P : 2 - P;
	
			if ( (i == 3 && i1 != 3) || (i1 == 3 && i != 3) )
				continue;

			x = (i < 3) ? -1.0 : 1.0;//исходное
			
			//x *= dx * dy * t_/ ( 16.0 * Pi * Pi);
			x *= dx * dy / ( 16.0 * Pi * Pi);
			
					
			for ( mf=c; mf <= M_fi; mf++ )
			{	
				//mr1 = mr;
				n = mf - 1;
									
				for ( mz=1; mz <= M_z; mz++ )
				{
					N_d = N_r[mz];
					n_hol = N_d-1;
					mb = m_b[N_d];
					mt = m_t[N_d];
						
					for ( mr=1; mr <= MR[n_hol]; mr++ )
					{
						nl = N_lin[n_hol][mr];//new
						m = (P) ? No2[i][mz][mr][mf] : No[i][mz][mr][mf];

						M = M_r_1;
						
						for ( mr1=1; mr1 <= M_r; mr1++ )
						{
							nl1 = N_lin[n_hol][mr1];//new
							if (abs(nl1 - nl) > 1)//new
								continue;
							
							for ( mz1=1; mz1 <= M_z; mz1++ )
							{
								if (N_d != N_r[mz1])
									continue;
								
								m1 = (P) ? No2[i1][mz1][mr1][mf] : No[i1][mz1][mr1][mf];

								if (abs(mz -mz1)> 1)
									continue;
								
								Rv = J[n][mr][n_hol];
								Rv1 = J[n][mr1][n_hol];
								
								Ry = (mr <= M) ? 0.0 : Y1[n][mr][n_hol];
								Ry1 = (mr1 <= M) ? 0.0 : Y1[n][mr1][n_hol];
								
								B = (nl == nl1) ? r_l[n_hol][nl] * Rv * Rv1 + r_l[n_hol][nl-1] * Ry * Ry1 :
										(nl == nl1 + 1) ? - r_l[n_hol][nl-1] * Ry * Rv1 : - r_l[n_hol][nl] * Rv * Ry1;
																								
								z = x * B * Ii(P,i,i1,n,mz);

								u = cos(2.0 * b_[mz] * ky) * cos(2.0 * b_[mz1] * ky);
								AL = Pi * (double)(mz-1)/ H_[N_d-1];
							
								if (mz == mz1)
									u2 = (mz == mt) ? hm[mz] / 3.0 
													: (mz == mb) ? hm[mz+1] / 3.0
																: (hm[mz] + hm[mz+1]) / 3.0;
								else
									if (mz == mz1 - 1)
										u2 = hm[mz+1] / 6.0;
									else
										if (mz == mz1 + 1)
											u2 = hm[mz] / 6.0;
										else
											u2 = 0.0;

								u = z * u2;

								u *= (mz > mb) ? (0.5 * Pi - atan( R_m / AL)) / AL : 1.0 / R_m;

								Alap[m][m1].real += u;
							}//mz1
						}//mr1
					}//mr
				}//mz
			}//mf
			
		}//i1
	}//i
}//CONVERG


double Ii(int P,int i,int i1,int m,int mz)//new
{
	int n,M,N_d;
	double h,S,Sm,C,Cm,f,R,pci,x,y,tet,sum,alf,bet,
		c,c_tet,s_tet,U,V,Z;

	M = 2 * m + 5;
	h = Pi / M;
	N_d = N_r[mz] - 1;

	for ( sum = 0.0, n=0; n < M; n++ )
	{
		pci = -0.5 * Pi + h * (0.5 + n);
		
		alf = C = cos(pci);
		bet = S = sin(pci);
				
		x = a_[mz] * (C + kx / R_m);
		y = b_[mz] * (S + ky / R_m);

		//x = a_[mz] * cos(pci);
		//y = b_[mz] * sin(pci);
				
		tet = atan2(y,x);

		//perehod k kosomu
		c_tet = cos(tet);
		s_tet = sin(tet);

		c = a_[mz] / b_[mz];

		U = c_tet * c_fi0[N_d] + c * s_tet * s_fi0[N_d];
		V = s_tet * c_fi0[N_d] - c_tet * s_fi0[N_d] /  c;
		
		tet = atan2(V,U);
		Z = sqrt(U * U + V * V);
		//perehod k kosomu
	
		tet -= FI_r_;
		
		R = 1.0 / (Z * sqrt(x * x + y * y));
		R *= R * R;

		f = tet * (double)m;
		Cm = cos(f);
		Sm = sin(f);
		
		
		
		f = tet * (double)m;
		if (P)
		{
			Cm = cos(f);
			Sm = sin(f);
		}
		else
		{
			Sm = cos(f);
			Cm = sin(f);
		}
		/*
		if ( i != i1)
			f = Sm * Cm * alf * bet;
		else
			f = (i == 1) ? Sm * Sm *alf * alf : (i == 2) ? Cm * Cm * bet * bet :  Cm * Cm;
		*/
		
		f = 0.0;
		if ( i != i1)
			f = Sm * Cm * alf * bet;
		else
			if (i!= 3)
				f = (i == 1) ? Sm * Sm * alf * alf : (i == 2) ? Cm * Cm * bet * bet :  Cm * Cm;
			
		
		sum += f * R;
	}//n

	return(sum * 4.0 * h / Pi);
}//Ii

