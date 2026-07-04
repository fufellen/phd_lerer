
#include "VIBRr.h"


void MATR_( int mz,int mz1)
{
	// i - index of observation component i=1,2,3 for Ex,Ey,Ex - corresp.
	// i1 - index of source component i1=1,2,3 for Ex,Ey,Ex - corresp.
	// c - start index for observation mf c=1 for cos(mf-1)Theta; c=2 for sin(mf-1)Theta [ see presentation p.6]
	// c1 - start index for source mf c1=1 for cos(mf-1)Theta; c1=2 for sin(mf-1)Theta [ see presentation p.6]
	// ii - index of result
	int c,c1,mf,mf1,mr,mr1,m,m1,bu,bu1,P,P1,i,i1,ii,Mr,Mr1,P_c,P_f,
		n_hol,n_hol1,mi,mi1,ie,ie1,mm_,mm1_,ms_,ms1_;
	double R,R1,rg,ig,u,v,u1,v1,r1,r2,ri1,ri2,gr,gi;
	MKL_Complex16 **ppA_Alap;
	
	ppA_Alap = Alap;
	
	n_hol = N_r[mz] - 1;
	Mr = MR[n_hol];
	n_hol1 = N_r[mz1] - 1;
	Mr1 = MR[n_hol1];
	
	P_c = 0;
	P_f = 1;

	for ( ii=0, i=1; i <= 3; i++ )
	{
		double ***AR,***AR1;
        int ***ANo,***ANo1;
      
		bu = (i == 1) ? 1 : 0;
		
		for ( i1=1; i1 <= 3; i1++ )
		{
			bu1 = (i1 == 1) ? 1 : 0;
			ii++;
		
			gr = G_r[ii];
			gi = G_i[ii];

			//double greenFuncNorm = gr * gr + gi * gi;
            if(gr * gr + gi * gi < 1e-5)
                continue;

			for ( P=0; P <= 1; P++ )
			{
				//double **g_r_z2 = g_r[ii];
               // double **g_i_z2 = g_i[ii];

				c = (bu) ? 1 + P: 2 - P;
				
				if (P)
				{
					AR = (bu)? R_s2 : R_c2;
						
					ANo = No2[i];//наблюдения observation
				
				}
				else
				{
					AR = (bu)? R_c2 : R_s2;
						
					ANo = No[i];//наблюдения observation
										
				}
				
				//MP_c = (i == i1) ? P : 1;//Proveril
				//MP_c = 1;
				
				//for ( P1=0; P1 <= MP_c; P1++ )
			//	if (BU_SIM)
			//		P_c = P_f = P;
				
				//P1 = P;
				//for ( P1=P_c; P1 <= P_f; P1++ )
				for ( P1=0; P1 <= 1; P1++ )
				{
					if (P1)
					{
						AR1 = (bu1)? R_s2 : R_c2;
						
						ANo1 = No2[i1];//наблюдения observation
					}
					else
					{
						AR1 = (bu1)? R_c2 : R_s2;
						
						ANo1 = No[i1];//наблюдения observation
										
					}

					c1 = (bu1) ? 1 + P1: 2 - P1;
						
					//for ( mz=1; mz <= ms; mz++ )
					{ // loop over observation
						//double *g_r_z = g_r_z2[mz];
                       // double *g_i_z = g_i_z2[mz];
						int **m_=ANo[mz];
						
						MKL_Complex16 *A_Alap;
								
						for ( mf=c; mf <= M_fi; mf++ )
						{
							
							double *AR_ = AR[mz][mf];

							for ( mr=1; mr <= Mr; mr++ )//new
							{
									
								R = AR_[mr];
								m = m_[mr][mf];
								A_Alap = ppA_Alap[m];
	
								//for ( mz1=1; mz1 <= ms; mz1++ )
								{
									int **m1_=ANo1[mz1];
										
									//rg = R * g_r_z[mz1];
									//ig = R * g_i_z[mz1];
									rg = R * gr;
									ig = R * gi;
																	
									for ( mf1=c1; mf1 <= M_fi; mf1++ )
									{
										double *AR1_ = AR1[mz1][mf1];

										for ( mr1=1; mr1 <=  Mr1; mr1++ )//new
										{	
											R1 = AR1_[mr1];
											m1 = m1_[mr1][mf1];
										
											u = R1 * rg;
											v = R1 * ig;
											A_Alap[m1].real += u;
											A_Alap[m1].imag += v;

											if (N_ell > 1)
											{
												for ( ie=2; ie <= N_ell; ie++ )//неоднородная часть The inhomogeneous part
												{
													//double *er_ = er[ie],*ei_ = ei[ie];
													mi = MS_* (ie - 1);
													
													mm_ = m + mi;
													ms1_ = m1 + mi;

													for (  ie1=1; ie1 < ie; ie1++ )
													{
														mi1 = MS_* (ie1 - 1);
														
														ms_ = m + mi1;
														mm1_ = m1 + mi1;
																											
														u1 = er[ie][ie1];
														v1 = ei[ie][ie1];
														//u1 = v1 = 0.0;//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
														r1 = u * u1;
														r2 = v * v1;
														ri1 = v * u1;
														ri2 = u * v1;
	
														Alap[mm_][mm1_].real += r1 - r2;
														Alap[mm_][mm1_].imag += ri1 + ri2;

														Alap[ms_][ms1_].real += r1 + r2;
														Alap[ms_][ms1_].imag += ri1 - ri2;

													}//ie1
													
												}//ie
											}//N_ell > 1
																				
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