
#include "VIBRr.h"


void PRO_SIM( void)
{
	// i - index of observation component i=1,2,3 for Ex,Ey,Ex - corresp.
	// i1 - index of source component i1=1,2,3 for Ex,Ey,Ex - corresp.
	// c - start index for observation mf c=1 for cos(mf-1)Theta; c=2 for sin(mf-1)Theta [ see presentation p.6]
	// c1 - start index for source mf c1=1 for cos(mf-1)Theta; c1=2 for sin(mf-1)Theta [ see presentation p.6]
	// ii - index of result
	int c,c1,mz,mz1,mf,mf1,mr,mr1,m,m1,bu,bu1,P,P1,i,i1,ii,MP_c,
	n_hol,n_hol1;
	//double R,R1,rg,ig;
//	MKL_Complex16 **ppA_Alap;
	
//	 ppA_Alap = Alap;
	
	
	for ( ii=0, i=1; i <= 3; i++ )
	{
		//double ***AR,***AR1;
        int ***ANo,***ANo1;
      
		bu = (i == 1) ? 1 : 0;
		
		for ( i1=1; i1 <= 3; i1++ )
		{
			bu1 = (i1 == 1) ? 1 : 0;
			ii++;
			
			for ( P=0; P <= 1; P++ )
			{
		//		double **g_r_z2 = g_r[ii];
         //       double **g_i_z2 = g_i[ii];

				c = (bu) ? 1 + P: 2 - P;
				
				if (P)
				{
				//	AR = (bu)? R_s2 : R_c2;
						
					ANo = No2[i];//наблюдения observation
				
				}
				else
				{
				//	AR = (bu)? R_c2 : R_s2;
						
					ANo = No[i];//наблюдения observation
										
				}
				
				//MP_c = (i == i1) ? P : 1;//Proveril
				MP_c = 1;
				
				for ( P1=0; P1 <= MP_c; P1++ )
				{
					if (P1)
					{
					//	AR1 = (bu1)? R_s2 : R_c2;
						
						ANo1 = No2[i1];//наблюдения observation
					}
					else
					{
					//	AR1 = (bu1)? R_c2 : R_s2;
						
						ANo1 = No[i1];//наблюдения observation
										
					}

					c1 = (bu1) ? 1 + P1: 2 - P1;
			
			
					for ( mz=1; mz <= ms; mz++ )
					{ // loop over observation
						//double *g_r_z = g_r_z2[mz];
                       // double *g_i_z = g_i_z2[mz];
						int **m_=ANo[mz];
						
					//	MKL_Complex16 *A_Alap;
						
						n_hol = N_r[mz] - 1;//new
							
						for ( mf=c; mf <= M_fi; mf++ )
						{
							
							//double *AR_ = AR[mz][mf];

							for ( mr=1; mr <= MR[n_hol]; mr++ )//new
							{
									
							//	R = AR_[mr];
								m = m_[mr][mf];
								// A_Alap = ppA_Alap[m];
	
								for ( mz1=1; mz1 <= ms; mz1++ )
								{
									int **m1_=ANo1[mz1];

									n_hol1 = N_r[mz1] - 1;//new
										
								//	rg = R * g_r_z[mz1];
								//	ig = R * g_i_z[mz1];
																		
									for ( mf1=c1; mf1 <= M_fi; mf1++ )
									{
									//	double *AR1_ = AR1[mz1][mf1];

										for ( mr1=1; mr1 <=  MR[n_hol1]; mr1++ )//new
										{	
									//		R1 = AR1_[mr1];
											m1 = m1_[mr1][mf1];
										
									//		A_Alap[m1].real += R1 * rg;
									//		A_Alap[m1].imag += R1 * ig;

											if (Alap[m][m1].real  ==Alap[m1][m].real && Alap[m][m1].imag  ==Alap[m1][m].imag )
											{
												printf ("\n\n P%3d   P1%3d    i%3d  i1%3d    mz%3d   mz1%3d    nh%3d  nh1%3d\n     mr%3d   mr1%3d  mf%3d   mf1%3d",
													P,P1,i,i1,mz,mz1,n_hol,n_hol1,mr,mr1,mf,mf1);   
		
												fprintf( out,"\n\n P%3d   P1%3d    i%3d  i1%3d    mz%3d   mz1%3d    nh%3d  nh1%3d\n     mr%3d   mr1%3d  mf%3d   mf1%3d",
													P,P1,i,i1,mz,mz1,n_hol,n_hol1,mr,mr1,mf,mf1);   
											}
																				
										}//mr1
									}//mf1
								}//mz1
								
							}//mr
						}//mf
								
					}//mz // loop over observation
				}//P1
			}//P
		}//i1
	}//i
			                                                                                                                                                             
}//MATR