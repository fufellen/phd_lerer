
#include "VIBRr.h"


//void MATR( int c,int c1,int i,int i1,int ii)
void MATR( int i,int i1,int ii)
{
	// i - index of observation component i=1,2,3 for Ex,Ey,Ex - corresp.
	// i1 - index of source component i1=1,2,3 for Ex,Ey,Ex - corresp.
	// c - start index for observation mf c=1 for cos(mf-1)Theta; c=2 for sin(mf-1)Theta [ see presentation p.6]
	// c1 - start index for source mf c1=1 for cos(mf-1)Theta; c1=2 for sin(mf-1)Theta [ see presentation p.6]
	// ii - index of result
	int c,c1,mz,mz1,mf,mf1,mr,mr1,m,m1,bu,bu1,P;
	double R,R1,rg,ig;
	
	bu = (i == 1) ? 1 : 0;
	bu1 = (i1 == 1) ? 1 : 0;
	
	for ( P=0; P <= 1; P++ )
	{
		//c = 2 - P;
		//c1 = (ii == 8) ? c : 1 + P;
		c = (ii <= 3) ? 1 + P: 2 - P;
		c1 = (ii == 1 || ii == 4 || ii == 7) ? 1 + P: 2 - P;

		for ( mr=1; mr <= M_r; mr++ )
			for ( mf=c; mf <= M_fi; mf++ )
			{ // loop over observation
				for ( mr1=1; mr1 <= M_r; mr1++ )
					for ( mf1=c1; mf1 <= M_fi; mf1++ )
					{ // loop over source
						for ( mz=1; mz <= ms; mz++ )
						{
							//R = (i == 1)? R_s[mz][mf][mr] : R_c[mz][mf][mr];
							if (P)
							{
								R = (bu)? R_s2[mz][mf][mr] : R_c2[mz][mf][mr];
						
								m = No2[i][mr][mf][mz];//наблюдения observation
							}
							else
							{
								R = (bu)? R_s[mz][mf][mr] : R_c[mz][mf][mr];
						
								m = No[i][mr][mf][mz];//наблюдения observation
							}
						
							//m = No[i][mr][mf][mz];//наблюдения observation

							for ( mz1=1; mz1 <= ms; mz1++ )
							{											
								rg = R * g_r[ii][mz][mz1];
								ig = R * g_i[ii][mz][mz1];
								if (P)
								{
									R1 = (bu1)? R_s2[mz1][mf1][mr1] : R_c2[mz1][mf1][mr1];
									
									m1 = No2[i1][mr1][mf1][mz1];//исток source

																		
									Alap2[m][m1].real += R1 * rg;
									Alap2[m][m1].imag += R1 * ig;
								}
								else
								{
									R1 = (bu1)? R_s[mz1][mf1][mr1] : R_c[mz1][mf1][mr1];
									
									m1 = No[i1][mr1][mf1][mz1];//исток source

																		
									Alap[m][m1].real += R1 * rg;
									Alap[m][m1].imag += R1 * ig;
								}
							}//mz1
						}//mz
					}//mf1 // loop over source
			}//mf // loop over observation
	}//P
			                                                                                                                                                             
}//MATR