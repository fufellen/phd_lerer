
#include "VIBRr.h"


void MATR_S( void)
{
	// i - index of observation component i=1,2,3 for Ex,Ey,Ex - corresp.
	// i1 - index of source component i1=1,2,3 for Ex,Ey,Ex - corresp.
	// c - start index for observation mf c=1 for cos(mf-1)Theta; c=2 for sin(mf-1)Theta [ see presentation p.6]
	// c1 - start index for source mf c1=1 for cos(mf-1)Theta; c1=2 for sin(mf-1)Theta [ see presentation p.6]
	// ii - index of result
	int c,c1,mz,mz1,mf,mf1,mr,mr1,m,m1,P,P1,i,i1,ii,m_,m1_,M_r_c;
		
	
	for ( ii=0, i=1; i <= 3; i++ )
	{
		for ( i1=1; i1 <= 3; i1++ )
		{
			ii++;

			if (i != i1)
				continue;
			for ( P=0; P <= 1; P++ )
			{
				c = (i == 1) ? 1 + P: 2 - P;
				for ( P1=0; P1 <= P; P1++ )
				{
					//P1 = P;
					//c1 = c;
				c1 = (i1 == 1) ? 1 + P1: 2 - P1;

				for ( mr=1; mr <= M_r; mr++ )
					{
						M_r_c = (i == i1 && P == P1) ? mr - 1 : M_r - 1;
					
						for ( mf=c; mf <= M_fi; mf++ )
						{ // loop over observation
					
							for ( mz=1; mz <= ms; mz++ )
							{								
								for ( mr1=1; mr1 <= M_r_c; mr1++ )
								for ( mf1=c1; mf1 <= M_fi; mf1++ )
											
								for ( mz1=1; mz1 <= ms; mz1++ )
								{											
									m = (P) ? No2[i][mr][mf][mz] : No[i][mr][mf][mz];
									m1 = (P1) ? No2[i1][mr1][mf1][mz1] : No[i1][mr1][mf1][mz1];
									
									m_ = (P1) ? No2[i][mr1][mf1][mz] : No[i][mr1][mf1][mz];
									m1_ = (P) ? No2[i1][mr][mf][mz1] : No[i1][mr][mf][mz1];
									
																					
									Alap[m_][m1_].real = Alap[m][m1].real;
									Alap[m_][m1_].imag = Alap[m][m1].imag;
					
								}//mz1
								
							}//mz
								
						}//mf // loop over observation
					}//mr
				}//P1
			}//P
		}//i1
	}//i
		
	
}//MATR