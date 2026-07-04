
#include "VIBRr.h"


void D_N_PER( void )//рассеянная мощность Scattered power
{
	int N_d,c1,i,i1,mz1,mz1_,mf1,mr1,q,p,m1,g,BU[3],q1,q2,p1,p2,Pe,Mrg,n_hol,ie,m1_,nl;//new!!!!!!!!
	double Rh,R1,x,y,alf,bet,gam2,Ro,eta,P[3],P0[3],Pg[3],alf2,bet2,//new!!!!!!!!
			u,v,ue,ve,uh,vh,***f,***fi,***fh,***fhi,del[3],Um1[3],Um2[3],
			//T_ex,
			u3,v3,
			***E_r_p,***H_r_p,***E_i_p,***H_i_p,P_p[3],Pg_p[3],sz[3],cz[3],
			G2i[111],G2[111],al,
			F[4],Fi[4],cgam[111],cgam_i[111],Cgam[111],Cgam_i[111],ga,ga2,ga_i,ga2_i,u1,u2,v1,v2,
			te[3],qu[3],fi1[3],fi3[3],S_1[3],S_3[3],GAM2[3][221],**S1,**C,
			u_,v_;
	
	
	//for ( Pol=0; Pol <= 1; Pol++ )
	{
	for ( Mrg=0, p1=0; p1 <= 1; p1++ )
		for ( p2=p1; p2 <= N_g; p2++ )
		{
			p = (p1) ? -p2 : p2;
		//	if (p)
		//		continue;
			
			alf = Pi_2 * (double)p / dx + kx;
			al = Pi_2 * (double)p / dx;
			alf2 = alf * alf;
			
			for ( q1=0; q1 <= 1; q1++ )
			for ( q2=q1; q2 <= N_g; q2++ )
			{
				q = (q1) ? -q2 : q2;
				//if (q)
					//continue;
			
				//bet = Pi_2 * (double)q / dy + ky - alf * tg_pci;
				bet = Pi_2 * (double)q / dy + ky - al * tg_pci;
				bet2 = bet * bet;
				Ro = alf2 + bet2;
				gam2 = Ro - k2_r[1];//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				
				BU[1] = (gam2 > 0.0) ? 1 : 0;
				
				// Ro - k2_r[5] > 0.0 this is a case of evanescent floquet harmonic (p,q) in lower halfspace k2_r[5]
				// otherwise Ro - k2_r[5] <= 0.0 is a case of propagating floquet harmonic (p,q) in lower halfspace k2_r[5]
				BU[2] = (Ro - k2_r[Nsl + 1] > 0.0 || EPS_i[Nsl + 1] < -1.0e-7 ) ? 1 : 0;

				//if (gam2 > 0.0)
				if (BU[1] && BU[2])
					continue;
					
				Mrg++;

			}//q
		}//p
	
		// integral of base function over (X,Y) of elliptic shape [see "presentation" p.6] 
	// accounting for symmetry of illumination that appears as cos or sin as R_c or R_s instead of exp depending on polarization


	f = ( double***)amat_3d( 1,2,1,3, 1,ms, sizeof( double) );
	fh = ( double***)amat_3d( 1,2,1,3, 1,ms, sizeof( double) );
	fi = ( double***)amat_3d( 1,2,1,3, 1,ms, sizeof( double) );
	fhi = ( double***)amat_3d( 1,2,1,3, 1,ms, sizeof( double) );
	S1 = ( double**)amat_2d( 1,2, 0,ms+1, sizeof( double) );
	C = ( double**)amat_2d( 1,2, 0,ms+1, sizeof( double) );

	E_r_p = ( double***)amat_3d( 1,Mrg,1,2,1,3,  sizeof( double) );
	E_i_p = ( double***)amat_3d( 1,Mrg,1,2,1,3,  sizeof( double) );
	H_r_p = ( double***)amat_3d( 1,Mrg,1,2,1,3,  sizeof( double) );
	H_i_p = ( double***)amat_3d( 1,Mrg,1,2,1,3,  sizeof( double) );

	for ( Pol=0; Pol <= 1; Pol++ )
	{
	for ( g=1; g <= Mrg; g++ )
		for ( i=1; i <= 2; i++ )
			for ( i1=1; i1 <= 2; i1++ )
				E_r_p[g][i1][i] = H_r_p[g][i1][i] = E_i_p[g][i1][i] = H_i_p[g][i1][i] = 0.0;
			
//	for ( Pe=0; Pe <= 1; Pe++ )
	{		
		// power density of downward field from source
		eta = k / Z0;
	
		for ( Mrg=0, p1=0; p1 <= 1; p1++ )
		for ( p2=p1; p2 <= N_g; p2++ )
		{
			p = (p1) ? -p2 : p2;
			
			alf = Pi_2 * (double)p / dx + kx;
			al = Pi_2 * (double)p / dx;
			alf2 = alf * alf;
			
			for ( q1=0; q1 <= 1; q1++ )
			for ( q2=q1; q2 <= N_g; q2++ )
			{
				q = (q1) ? -q2 : q2;
				//bet = Pi_2 * (double)q / dy + ky - alf * tg_pci;
				bet = Pi_2 * (double)q / dy + ky - al * tg_pci;
				
				bet2 = bet * bet;
				Ro = alf2 + bet2;
				gam2 = Ro - k2_r[1];//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

				
				for ( nl=1; nl <= N_ell; nl++ )
				{
					u = alf * XC[nl] + bet * YC[nl];
					Co_c[nl] = cos(u);
					Si_c[nl] = sin(u);
				}//nl

				// gam2 > 0.0 this is a case of evanescent floquet harmonic (p,q) in upper halfspace k2_r[1]
				// otherwise gam2 <= 0.0 is a case of propagating floquet harmonic (p,q) in upper halfspace k2_r[1]
				BU[1] = (gam2 > 0.0) ? 1 : 0;


				// Ro - k2_r[5] > 0.0 this is a case of evanescent floquet harmonic (p,q) in lower halfspace k2_r[5]
				// otherwise Ro - k2_r[5] <= 0.0 is a case of propagating floquet harmonic (p,q) in lower halfspace k2_r[5]
				//BU[2] = (Ro - k2_r[Nsl + 1] > 0.0) ? 1 : 0;
				BU[2] = (Ro - k2_r[Nsl + 1] > 0.0 || EPS_i[Nsl + 1] < -1.0e-7 ) ? 1 : 0;

				if (BU[1] && BU[2])
					continue;
				
				BASIS(alf,bet);//вычисление Rnp calculate//new
					
			//	printf ("\n  1  %1d  %1d ",p,q);                                                      
			//	fprintf( out,"\n  1  %1d  %1d ",p,q);   

				Mrg++;

				for ( i=1; i <= Nsl + 1; i++ )
				{
					u = G2[i] = Ro - k2_r[i];
					v = G2i[i] =  -k2_i[i];
					sq_compl(u,v, &u1,&v1);
				
					Cgam[i] = u1;//gamma
					Cgam_i[i] = v1;
		
					cgam[i] = v1; //gamma1 = i * gamma
					cgam_i[i] = -u1;
				}

				for ( mz1=1; mz1 <= ms; mz1++ )
				{
					N_d = N_r[mz1];

					GAM2[1][N_d] = Ro - k2_r[N_d];
					GAM2[2][N_d] = -k2_i[N_d];

					x = (Zm[mz1] - h_b[N_d]) * Cgam[N_d];
					y = (Zm[mz1] - h_b[N_d]) * Cgam_i[N_d];
				
					sinh_C(x,y, &u,&v);
					S1[1][mz1] = u;
					S1[2][mz1] = v;

					cosh_C(x,y, &u,&v);
					C[1][mz1] = u;
					C[2][mz1] = v;
				}//mz1
				//____________________________________________________________________________________
			
									
				/*
					E_[1][i] - отраженная волна reflected wave
					E_[2][i] - прошедшая волна transmitted wave
				*/

				for ( g=1; g <= 2; g++ )
				{ // g=1 for reflected field; g=2 for propagating field
					if (BU[g])
						continue;

					if ( g == 1)
					{
						ga = cgam[1]; 
						ga_i = cgam_i[1];
						ga2 = cgam[1] / EPS_r[1]; 
						ga2_i = cgam_i[1] / EPS_r[1];;
					}
					else
					{
						ga = -cgam[Nsl + 1]; 
						ga_i = -cgam_i[Nsl + 1];
				
						divide_(ga,ga_i, EPS_r[Nsl + 1],EPS_i[Nsl + 1],&u,&v );
						ga2 = u;
						ga2_i = v;
					}				

					for ( mz1=1; mz1 <= ms; mz1++ )//W
					{
						//mz1_ = ms - (mz1 - 1);ISPRAV

						N_d = N_r[mz1];
						mz1_ = m_t[N_d] - (mz1 - m_b[N_d]);
				//		if (N_d != N_r[mz1_])
				//			printf ("\n OSHIBKA 8\n");
										
						for ( i=1; i <= 3; i++ )
						{
							CONST(N_d,i,Ro, te,qu,fi1,fi3,del,S_1,S_3);
											

							x = Cgam[N_d] * H_[N_d-1];
							y = Cgam_i[N_d] * H_[N_d-1];			
							sinh_C(x,y, &u1,&v1);// sh(gamma * h)
							divide_(del[1],del[2], u1,v1, &u,&v);// delta_1 /sh(gamma * h)
						
							INT(N_d,1,mz1, Cgam[N_d], Cgam_i[N_d], S1,C, sz,cz);
							mult_(sz[1],sz[2],u,v, &Um1[1],&Um1[2]);//Um1
						
							INT(N_d,1,mz1_, Cgam[N_d], Cgam_i[N_d], S1,C, sz,cz);
							mult_(sz[1],sz[2],u,v, &Um2[1],&Um2[2]);//Um2
						
							if (g == 1)
							{
								mult_(Um1[1],Um1[2],fi3[1],fi3[2],&u,&v);
								mult_(Um2[1],Um2[2],qu[1],qu[2] ,&u1,&v1);
								u += u1;
								v += v1;

								mult_(u,v, S_1[1],S_1[2],&u,&v);
							}
							else
							{
								mult_(Um2[1],Um2[2],fi1[1],fi1[2],&u,&v);
								mult_(Um1[1],Um1[2],qu[1],qu[2] ,&u1,&v1);
								u += u1;
								v += v1;

								mult_(u,v, S_3[1],S_3[2],&u,&v);
								u = - u;
								v = - v;
							}
						
							u1 = u;
							v1 = v;//
						
							if (i != 3)
							{
								F[i] = u1 / Ro;
								Fi[i] = v1 / Ro;
							}
							else
							{
								mult_(u1,v1,EPS_r[N_d],EPS_i[N_d], &u1,&v1);
								F[i] = u1;
								Fi[i] = v1;
							}
						
						}//i
	//E-field
						//alf = - alf;
						mult_(F[1],Fi[1], k2_r[N_d],k2_i[N_d], &u1,&v1);//F1 c shlap
						mult_(F[2],Fi[2], GAM2[1][N_d],GAM2[2][N_d], &u2,&v2);//F2 c shlap
					
						f[1][1][mz1] = bet2 * u1 - alf2 * u2;
						fi[1][1][mz1] = bet2 * v1 - alf2 * v2;

						
						f[1][2][mz1] = f[2][1][mz1] = -alf * bet * (u1 + u2);//!!!!!!
						fi[1][2][mz1] = fi[2][1][mz1] = -alf * bet * (v1 + v2);//!!!!!!
						
						f[2][2][mz1] = alf2 * u1 - bet2 * u2;
						fi[2][2][mz1] = alf2 * v1 - bet2 * v2;
	
						mult_(ga2,ga2_i, F[3],Fi[3], &u3,&v3); 

						f[1][3][mz1] = alf * v3;
						fi[1][3][mz1] = -alf * u3;

						f[2][3][mz1] = bet * v3;
						fi[2][3][mz1] = -bet * u3;

						//alf = - alf;
	//END E-field

	//H-field
						//alf = - alf;
						mult_(F[1],Fi[1], EPS_r[N_d],EPS_i[N_d], &u1,&v1);//eps2 F1=F1 s chertoi
			
						if (g == 1)
							divide_(EPS_r[1],EPS_i[1],G2[1],G2i[1], &u2,&v2);
						else
							divide_(EPS_r[Nsl + 1],EPS_i[Nsl + 1],G2[Nsl + 1],G2i[Nsl + 1], &u2,&v2);
					
						mult_(G2[N_d],G2i[N_d], u2,v2, &u2,&v2);
						mult_(F[2],Fi[2], u2,v2, &u2,&v2);//F2 = F2 * gam2^2 * epsi / gami^2 = 
					
						mult_(ga,ga_i, u2 - u1,v2 - v1, &u,&v);
						fh[1][1][mz1] = -alf * bet * u;//!!!!!!
						fhi[1][1][mz1] = -alf * bet * v;//!!!!!!
					
						mult_(ga,ga_i, bet2 * u2 + alf2 * u1,bet2 * v2 + alf2 * v1, &u,&v);
						fh[1][2][mz1] = -u;
						fhi[1][2][mz1] = -v;
										
						mult_(ga,ga_i, alf2 * u2 + bet2 * u1,alf2 * v2 + bet2 * v1, &u,&v);
						fh[2][1][mz1] = u;
						fhi[2][1][mz1] = v;
			
						fh[2][2][mz1] = -fh[1][1][mz1];
						fhi[2][2][mz1] = -fhi[1][1][mz1];

						fh[1][3][mz1] = -bet * Fi[3];
						fhi[1][3][mz1] = bet * F[3];
					
						fh[2][3][mz1] = alf * Fi[3];
						fhi[2][3][mz1] = -alf * F[3];
				
						//alf = - alf;
	//END H-field
						
//_____________________________________________________________
						
					}//mz1 end W
					
					for ( Pe=0; Pe <= 1; Pe++ )
					for ( i=1; i <= 2; i++ )
					{
						for ( i1=1; i1 <= 3; i1++ )
						{
							c1 = (i1 == 1) ? 1 + Pe : 2 - Pe;
							//new							
							for ( mf1=c1; mf1 <= M_fi; mf1++ )
										
								for ( mz1=1; mz1 <= ms; mz1++ )
								{
									n_hol = N_r[mz1]-1;
								
									for ( mr1=1; mr1 <= MR[n_hol]; mr1++ )
									{
										//new
										if (Pe)
											R1 = (i1 == 1)? R_s2[mz1][mf1][mr1] : R_c2[mz1][mf1][mr1];
										else
											R1 = (i1 == 1)? R_c2[mz1][mf1][mr1] : R_s2[mz1][mf1][mr1];//new
											//R1 = (i1 == 1)? R_s[mz1][mf1][mr1] : R_c[mz1][mf1][mr1];
									
										//new end
									
										R1 *= 0.5;
										Rh = R1;
									
	
										m1 = (Pe) ? No2[i1][mz1][mr1][mf1] : No[i1][mz1][mr1][mf1];//исток
								
										u = f[i][i1][mz1];
										v = fi[i][i1][mz1];

										u_ = fh[i][i1][mz1] * eta;
										v_ = fhi[i][i1][mz1] * eta;

										for ( ie=1; ie <= N_ell; ie++ )
										{
											m1_ = m1 + MS_* (ie - 1);
											/*
											u1 = ASR_[Pol][m1+1];
											v1 = ASI_[Pol][m1+1];
											*/
											u1 = ASR_[Pol][m1_+1];
											v1 = ASI_[Pol][m1_+1];

											mult_(u1,v1, Co_c[ie], -Si_c[ie], &u1,&v1);

											mult_(-u,-v, u1,v1, &ue,&ve);//P * X
								
											//u = fh[i][i1][mz1] * eta;
											//v = fhi[i][i1][mz1] * eta;
											mult_(-u_,-v_, u1,v1, &uh,&vh);//P * X
								
											if (g == 1)
											{
												vh = - vh;
												ve = - ve;
												uh = - uh;
												ue = - ue;
											}
									
											E_r_p[Mrg][g][i] += R1 * ue; //дифрагированное поле diffracted field
											E_i_p[Mrg][g][i] += R1 * ve;

											H_r_p[Mrg][g][i] += Rh * uh; 
											H_i_p[Mrg][g][i] += Rh * vh;
										}//ie
									}//mr1
						
								}//mz1
							
							
						}//i1
					}//i
				}//g

				/*
				E[N_d][n][i][j], H[N_d][n][i][j]
				n = 1 -действительная часть, n = 2 - мнимая  частью
				i = 1 - x -овая компонента Е (H), i = 2 - e -овая компонента Е(H). 
				j = 1 - коэффициент отраженния соответствующей компоненты поля,
				j = 2 - коэффициент прохождения соответствующей компоненты поля,
				j = 3 - поле внутреннее с коэффициентом А,
				j = 4 - поле внутреннее с коэффициентом B.
				n = 1 the real part, n = 2 - imaginary part
				 i = 1 - x-EW component E (H), i = 2 - e-marketing component E (H).
				 j = 1 - reflectance of the corresponding component of the field,
				 j = 2 - transmission coefficient of the corresponding components of the field,
				 j = 3 - internal field with a coefficient A,
				 j = 4 - the internal field by a factor of B.
				 							
				*/
			}//q
		}//p
	}//Pe	

	N_d = 2;

	P0[0] = P[0] = Ey_0 * Hx_0 - Ex_0 * Hy_0;
	{//******
		double Er[3][4],Ei[3][4],Hr[3][4],Hi[3][4],Er_[3][4],Ei_[3][4],Hr_[3][4],Hi_[3][4],Pg_p_p[3],Pg_p_s[3];
		
		P_p[1] = P_p[2] = P[1] = P[2] = P0[1] = P0[2] = 0.0;
	
		for ( Mrg=0, p1=0; p1 <= 1; p1++ )
		for ( p2=p1; p2 <= N_g; p2++ )
		{
			p = (p1) ? -p2 : p2;
			
			alf = Pi_2 * (double)p / dx + kx;
			al = Pi_2 * (double)p / dx;
			alf2 = alf * alf;
			
			for ( q1=0; q1 <= 1; q1++ )
			for ( q2=q1; q2 <= N_g; q2++ )
			{
				q = (q1) ? -q2 : q2;
				bet = Pi_2 * (double)q / dy + ky - al * tg_pci;
				bet2 = bet * bet;
				Ro = alf2 + bet2;
				gam2 = Ro - k2_r[1];//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				
				BU[1] = (gam2 > 0.0) ? 1 : 0;
				
				// Ro - k2_r[5] > 0.0 this is a case of evanescent floquet harmonic (p,q) in lower halfspace k2_r[5]
				// otherwise Ro - k2_r[5] <= 0.0 is a case of propagating floquet harmonic (p,q) in lower halfspace k2_r[5]
				BU[2] = (Ro - k2_r[Nsl + 1] > 0.0 || EPS_i[Nsl + 1] < -1.0e-7 ) ? 1 : 0;
				
				if (BU[1] && BU[2])
					continue;
				
			//	printf ("\n       2  %1d  %1d ",p,q);                                                      
			//	fprintf( out,"\n      2  %1d  %1d ",p,q);   

				Mrg++;	

				
				for ( g=1; g <= 2; g++ )
				{
					if (BU[g])
					{
						Pg_p[g] = Pg[g] = 0.0;
						continue;
					}
					for ( i=1; i <= 2; i++ )
					{
						Er[g][i] = E_r_p[Mrg][g][i];
						Ei[g][i] = E_i_p[Mrg][g][i];
						Hr[g][i] = H_r_p[Mrg][g][i];
						Hi[g][i] = H_i_p[Mrg][g][i];
					}//i

					//for ( i=1; i <= 2; i++ )
					//	Er[g][i] = Ei[g][i] = Hr[g][i] = Hi[g][i] = 0.0;//Проверка внешнего поля Правильно для 2-х поляризаций
						
					if ( !p && !q )
					{				
						
						if (Pol)
						{
							for ( i=1; i <= 2; i++ )
							{					
								Er[g][i] += Ep[N_d][1][i][g];//Ei-ref полное поле с учетом внешнего full field of view of external

								Ei[g][i] += Ep[N_d][2][i][g];
					
								Hr[g][i] += Hp[N_d][1][i][g];//Hi-refl
								Hi[g][i] += Hp[N_d][2][i][g];
							
							}//i
						}//Pol
						else
						{
							for ( i=1; i <= 2; i++ )
							{					
								
								Er[g][i] += E[N_d][1][i][g];//Ei-ref полное поле с учетом внешнего full field of view of external

								Ei[g][i] += E[N_d][2][i][g];
					
								Hr[g][i] += H[N_d][1][i][g];//Hi-refl
								Hi[g][i] += H[N_d][2][i][g];
							}//i
						}//!Pol
				
					}//!p && !q

					/**/
					mult_(Er[g][1],Ei[g][1], Hr[g][2],-Hi[g][2], &u,&v);// полная мощность гармоники full power harmonics

					Pg_p[g] = u / P[0];

					mult_(Er[g][2],Ei[g][2], Hr[g][1],-Hi[g][1], &u,&v);// полная мощность гармоники full power harmonics

					Pg_p[g] -= u / P[0];

					//Power(g,Er,Ei, Hr,Hi, Pg_p);
					//Pg_p[g] /= P[0];
					if ( !p && !q )
					{
						Er_[g][1] = Er[g][1] * cofi + Er[g][2] * sifi;//в плоскости падения 
						Ei_[g][1] = Ei[g][1] * cofi + Ei[g][2] * sifi;

						Er_[g][2] = Er[g][2] * cofi - Er[g][1] * sifi;//перпендикулярн плоскости падения
						Ei_[g][2] = Ei[g][2] * cofi - Ei[g][1] * sifi;

						Hr_[g][1] = Hr[g][1] * cofi + Hr[g][2] * sifi;//в плоскости падения
						Hi_[g][1] = Hi[g][1] * cofi + Hi[g][2] * sifi;
					
						Hr_[g][2] = Hr[g][2] * cofi - Hr[g][1] * sifi;//перпендикулярн плоскости падения
						Hi_[g][2] = Hi[g][2] * cofi - Hi[g][1] * sifi;

						mult_(Er_[g][1],Ei_[g][1], Hr_[g][2],-Hi_[g][2], &u,&v);//p - POLARIZATION
						Pg_p_p[g] = u / P[0];

						mult_(Er_[g][2],Ei_[g][2], Hr_[g][1],-Hi_[g][1], &u,&v);//s - POLARIZATION
						Pg_p_s[g] = -u / P[0];
					}//!p && !q
					
					/*
					mult_(E_r[g][1],E_i[g][1], H_r[g][2],-H_i[g][2], &u,&v);// только дифрагированная мощность гармоники only the diffracted power harmonics

					Pg[g] = u / P[0];

					mult_(E_r[g][2],E_i[g][2], H_r[g][1],-H_i[g][1], &u,&v);// только дифрагированная мощность гармоники only the diffracted power harmonics

					Pg[g] -= u / P[0];
					*/

					

					if (g == 2)
					
						Pg_p[2] = - Pg_p[2];
						//Pg[2] = - Pg[2];
					
					//if (p)
					//{
					//	Pg_p[g] /= 2.0;
					////	Pg[g] /= 2.0;
					//}
				
					P_p[g] += Pg_p[g];// полная мощность всех гармоник full power of all harmonics

				//	P[g] += Pg[g];//полная мощность только дифрагированная мощность всех гармоник full power only the diffracted power of all harmonics

				
					
				}//g

				if (!p &&!q)
				{
					if (Pol)
					{
						printf ("\n\n\n         p - POLARIZATION");   
						fprintf( out,"\n\n\n         p - POLARIZATION"); 
					}
					else
					{
						printf ("\n\n         s - POLARIZATION");   
						fprintf( out,"\n\n         s - POLARIZATION"); 
					}
			
					printf ("\n  %1d  %1d %7.5lf  %7.5lf",p,q ,Pg_p[1] ,Pg_p[2]  );                                                      
					fprintf( out,"\n  %1d  %1d %7.5lf  %7.5lf",p,q ,Pg_p[1] ,Pg_p[2]  );   
					if (graf2)
					{
						if (!Pol)
							fprintf(out2, "\n  %7.5lf    ", lam );
						fprintf(out2, "  %7.5lf  %7.5lf", Pg_p[1], Pg_p[2]);
					}


					//printf ("       %5.3e  ", P_p[1]);   
					//fprintf( out,"       %5.3e  ", P_p[1]); 
					//printf ("\n       %1d  %1d %5.3e  %5.3e",p,q ,sqrt(Pg[1]) ,sqrt(Pg[2])  );                                                      
					//fprintf( out,"\n      %1d  %1d %5.3e  %5.3e",p,q ,sqrt(Pg[1]) ,sqrt(Pg[2])  );
				}
				
			}//q
		}//p

		printf ("  %5.3e", 1.0 - (P_p[1] + P_p[2]) );
		fprintf( out," %5.3e",1.0 - (P_p[1] + P_p[2]) );
		if (graf2)
			fprintf(out2, " %5.3e", 1.0 - (P_p[1] + P_p[2]));

		printf ("\n		            	 To p: %7.5lf      To s: %7.5lf", Pg_p_p[1], Pg_p_s[1] );
		fprintf( out,"\n	            		 To p: %7.5lf      To s: %7.5lf", Pg_p_p[1], Pg_p_s[1] );

		

/*
	
		{
			T_ex = 1.0;
	
			printf ("\n           %5.3e  %5.3e      %5.3e  %5.3e", P_p[1], P_p[2], P[1], P[2]);   
			fprintf( out,"\n            %5.3e  %5.3e      %5.3e  %5.3e", P_p[1], P_p[2], P[1], P[2]);  
	}
*/		
		//printf ("   %5.3e", P[2] * T_ex );
	
		/*if (fabs(1-(P0[1] + P0[2])) > 1.0e-12) 
		{
			printf (" //%5.3e  %5.3e",sqrt(Pg0[1]) ,sqrt(Pg0[2])  );                                                      
			//fprintf( out," //%5.3e  %5.3e",sqrt(Pg0[1]) ,sqrt(Pg0[2])  );
			printf ("  %5.3e",1-(P0[1] + P0[2]) );
			//fprintf( out,"  %5.3e",1-(P0[1] + P0[2]) );
		}
		*/
	
		
	}//*****
	}//new Pol
	fmat_3d( (void***)fhi, 1,2, 1,3, 1, sizeof( double) );
	fmat_3d( (void***)fi, 1,2, 1,3, 1, sizeof( double) );
	fmat_3d( (void***)fh, 1,2, 1,3, 1, sizeof( double) );
	fmat_3d( (void***)f, 1,2, 1,3, 1, sizeof( double) );
	fmat_2d( (void**)S1, 1,2, 0, sizeof( double) );
	fmat_2d( (void**)C, 1,2, 0, sizeof( double) );

	fmat_3d( (void***)E_r_p, 1,Mrg, 1,2, 1, sizeof( double) );
	fmat_3d( (void***)E_i_p, 1,Mrg, 1,2, 1, sizeof( double) );
	fmat_3d( (void***)H_r_p, 1,Mrg, 1,2, 1, sizeof( double) );
	fmat_3d( (void***)H_i_p, 1,Mrg, 1,2, 1, sizeof( double) );
	
	}//Pol
	
}/* D_N_PER */

