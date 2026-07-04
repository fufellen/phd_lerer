
#include "VIBRr.h"


void F_df( int qu)// 
{ // !!!core function 
	 
	int i,i1,mz,mz1,q,q1,q2,p,m,n,buv,bur,mz_,mz1_,bur1,bur2,N_d1,N_d,p1,p2,
		m1,m_,m1_,p_g,q_g;
	double x,y,alf,alf2,bet,bet2,gam2,gam2_i,R_m,R_m2,u1,v1,u2,v2,ftol = 5.0,
			u,v,f[4][4],f_i[4][4],ar,br,abr,k2_2,k2_2_i,VE[3][6],U_c[3][3],
			bet1,**rq,**iq,ro_,***A_,***B_,D_1[3],D_[3],C_[3],s1[3],t1[3],
			M_1,M_11,SR=1.0,ft_e,
			fi1[3][4][221],fi3[3][4][221],Q[3][4][221],f1_[3],f3_[3],
			**rte,**ite,rg[4][111],ig[4][111],G2[3][111],
			del[4][221],deli[4][221],rga3[4][221],iga3[4][221],rga1[4][221],iga1[4][221],
			RM2,Rex=6.0,Ro,al,
			**rq_,**iq_,s2[3];
	//double um,vm,wm;

	if( kbhit() )
	   if( getch() == ESCape )
		 {
		 fflush( out );
		 exit( 1 );
		 }

	AM6();
	
	ft_e = - 5.0 *log(ft_l);

	rte = ( double**)amat_2d( 1,3, 1,Nsl+1, sizeof( double) );
	ite = ( double**)amat_2d( 1,3, 1,Nsl+1, sizeof( double) );
	rq = ( double**)amat_2d( 1,3, 1,Nsl+1, sizeof( double) );
	iq = ( double**)amat_2d( 1,3, 1,Nsl+1, sizeof( double) );
	rq_ = ( double**)amat_2d( 1,3, 1,Nsl+1, sizeof( double) );
	iq_ = ( double**)amat_2d( 1,3, 1,Nsl+1, sizeof( double) );
	
	A_ = ( double***)amat_3d( 1,2,1,3, 1,M_z, sizeof( double) );
	B_ = ( double***)amat_3d( 1,2,1,3, 1,M_z, sizeof( double) );

	RM2 = Rad_ro;
	p_g = (int)(dx / lam);
	q_g = (int)(dy / lam);
	for ( m=1; m <= 2; m++ )
	for ( mz=0; mz <= ms+1; mz++ )//наблюдения observation
	{
		Sh[m][mz] = Ch[m][mz] = Vd[m][mz] = 0.0;
		for ( mz1=0; mz1 <= ms+1; mz1++ )//исток  source
			U[m][mz][mz1] = V[m][mz][mz1] = Wu[m][mz][mz1] = Wv[m][mz][mz1] = 0.0;
	}
	
	N_corn = 2;
	for ( m=0; m < MS_N; m++ )
	{		
		for ( n=0; n < N_corn; n++ )
			Blap[n][m].real = Blap[n][m].imag = 0.0;

		for ( n=0; n < MS_N; n++ )
			Alap[m][n].real = Alap[m][n].imag = 0.0;
	}//m

	//_______________________________________________________________________
	NOMER( ); // maping of 4 indices[(X,Y,Z)l,m,n] of base function in eq.(2) of presentation into 1 index that becomes index of unknown in 
	
		
	//___________________________________________________________________
	////R_m = Max(Pi_2 * (double)N_g / dx,Pi_2 * (double)N_g / dy);//Huje pri dx!=dy!
	R_m = (Pi_2 * (double)N_g / dx + kx < Pi_2 * (double)N_g / dy + ky) ? Pi_2 * (double)N_g / dx  + kx 
																										: Pi_2 * (double)N_g / dy + ky;

	
	R_m2 = 0.5;
	R_m2 *= R_m;

	R_m *= R_m;//!!!!!!!!!!!!!!!!!
	R_m2 *= R_m2;

	DROITE();//// calc. Integral{Ee*V(mfi,mr,mz)}dV for external field Ee multiplied by Base function over elliptic cylinder volume.
	// result is in Blap i.e. right side of the equation system;
	
	for ( p1=0; p1 <= 1; p1++ )
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
			bet1 = bet;
			
			bet2 = bet * bet;
		
			Ro = alf2 + bet2;
			ro_ = sqrt(Ro);
	
			//if ( abs(p) > p_g  && abs(q) > q_g && Ro > R_m)
			if ( p  && q && Ro > R_m)
			//if ( Ro > R_m)
				break; // Ro > R_m case has been already treated in CONVERG
			
			for ( i=1; i <= N_ell; i++ )
			{
				u = alf * XC[i] + bet * YC[i];
				for ( i1=1; i1 < i; i1++ )
				{
					u1 = u - (alf * XC[i1] + bet * YC[i1]);
					er[i][i1] = cos(u1);
					ei[i][i1] = sin(u1);
				}
			}

			BASIS(alf,bet1);//вычисление Rnp calculate

			// further case Ro <= R_m
			buv = ( (!p && !q) || Ro <= R_m / RM2) ? 1 : 0;


			ar = alf2 / Ro;
			br = bet2 / Ro;
			abr = alf * bet / Ro;
						
			for ( n=1; n <= Nsl+1; n++ )
			{ 
				u2 = Ro - k2_r[n]; // GammaSquared in each layer

				if (bu_p[n])
				{
					u1 = -k2_i[n];
					sq_compl(u2,u1, &u,&v);
										
					divide_(1.0, 0.0, u2,u1, &Gm1[1][n],&Gm1[2][n]);
				}//bu_p[n]
				else
				{
					u1 = 0.0;
					Gm1[1][n] = 1.0 / u2;
					Gm1[2][n] = 0.0;

					if (u2 > 0.0)
					{
						u = sqrt(u2);
						v = 0.0;
					}
					else
					{
						v = sqrt(-u2);
						u = 0.0;
					}
				}//!bu_p[n]
	
				rg_[n] = u;//gamma
				ig_[n] = v;
					
			}//n

			if (!buv)
				goto met; 

			
//________________________________________________________________________________
			 // half-space upproximation for all underlying layers
			for ( n=1; n <= Nsl+1; n++ )//Вачисление gamma[1:3] с черточкой и q, t Vachislenie gamma [1:3] with a dash and q, t
			// see "Matrichnie elementi" p.3
			{ // loop over Nsl+1 layers
				G2[1][n] = Ro - k2_r[n]; // GammaSquared in each layer
				G2[2][n] = u1 = -k2_i[n];
					
				if (n > 1 && n <= Nsl) 
				{	
					u = rg_[n] * H_[n-1];
					v = ig_[n] * H_[n-1];
					SH_CTH_2(ro_ * H_[n-1],u,v,&s1[1],&s1[2],&s2[1],&s2[2],&t1[1],&t1[2]);
				}

				for ( i=1; i <= 3; i++ )
				{					
					if ( i == 1)
					{
						dze[i][n] = 1.0;
						dze_im[i][n] = 0.0;
					}
					else
						{
							if (i == 3)
								divide_(1.0,0.0, EPS_r[n],EPS_i[n], &x,&y);
							else
								divide_(EPS_r[n],EPS_i[n], G2[1][n],G2[2][n], &x,&y);
														
							dze[i][n] = x;
							dze_im[i][n] = y;
						}

					mult_(rg_[n],ig_[n], dze[i][n],dze_im[i][n], &u,&v);
					
					rg[i][n] = u;//gamma[1:3] с черточкой. Первый индекс номер компоненты ФГ.with a dash. The first index number of components of the FG
					ig[i][n] = v;

					if (n > 1 && n <= Nsl) 
					{	
						mult_(rg[i][n],ig[i][n], s1[1],s1[2], &u,&v);
						rq[i][n] = u;
						iq[i][n] = v;

						mult_(rg[i][n],ig[i][n], s2[1],s2[2], &u,&v);
						rq_[i][n] = u;
						iq_[i][n] = v;

						mult_(rg[i][n],ig[i][n], t1[1],t1[2], &u,&v);
						rte[i][n] = u;
						ite[i][n] = v;
												
					}//n > 1 &&  n < Nsl+1
					else
					{
						rq[i][n] = iq[i][n] = rq_[i][n] = iq_[i][n] = 0.0;

						rte[i][n] = rg[i][n];
						ite[i][n] = ig[i][n];
					}
				
				}//i
					
			}//n

//________________________________________________________________________________ !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                 
			
			for ( i=1; i <= 3; i++ )
				for ( N_d=2; N_d <= N_hol + 1; N_d++ )
				{
				//	if (abs(p)<2  && abs(q)<2)
						RECURS_ZERO(Ro,N_d,i,rte,ite,rq,iq, &u1,&v1, &u2,&v2);
				//	else
				//		RECURS(Ro,N_d,i,rte,ite,rq_,iq_, &u1,&v1, &u2,&v2);
					
					rga1[i][N_d] = u1;
					iga1[i][N_d] = v1;
					rga3[i][N_d] = u2;
					iga3[i][N_d] = v2;
				}
						
//________________________________________________________________________________
			met: 
			
			bur = 0;
						
			bur_v = (!bur || buv) ? 1 : 0;//ispolzuetsya v MATR
						
//_________________________________________
										
			Ancillary();

			if (!buv)
				goto met2;
			//____________________________________________________//calculate Вачисление dzeta / delta
			 // half-space upproximation for all underlying layers
			// see "Matrichnie elementi" p.2
			for ( i=1; i <= 3; i++ )
				for ( N_d=2; N_d <= N_hol + 1; N_d++ )
				{
					f1_[1] = rte[i][N_d] + rga1[i][N_d];
					f1_[2] = ite[i][N_d] + iga1[i][N_d];
					f3_[1] = rte[i][N_d] + rga3[i][N_d];
					f3_[2] = ite[i][N_d] + iga3[i][N_d];

					mult_(f1_[1],f1_[2], f3_[1],f3_[2], &x,&y);
					x -= rq[i][N_d] * rq[i][N_d] - iq[i][N_d] * iq[i][N_d];
					y -= 2.0 * rq[i][N_d] * iq[i][N_d];
				
					divide_(dze[i][N_d],dze_im[i][N_d],x,y, &del[i][N_d],&deli[i][N_d]); 
					mult_(f1_[1],f1_[2], del[i][N_d],deli[i][N_d], &fi1[1][i][N_d],&fi1[2][i][N_d]);
					mult_(f3_[1],f3_[2], del[i][N_d],deli[i][N_d], &fi3[1][i][N_d],&fi3[2][i][N_d]);
					mult_(rq[i][N_d],iq[i][N_d], del[i][N_d],deli[i][N_d], &Q[1][i][N_d],&Q[2][i][N_d]);
					
				}//N_d, i
		
			//____________________________________________________
			met2: 
					
			for ( i=1; i <= 3; i++ )//неоднородная часть The inhomogeneous part
				for ( i1=1; i1 <= 3; i1++ )
					f_i[i][i1] = 0.0;
		
			f[1][2] = f[2][1] = -alf * bet;
			f[1][3] = alf;
			f[2][3] = bet;
			f[3][1] = -f[1][3];
			f[3][2] = -f[2][3];
			f[3][3] = Ro;
			

			//____________________________________________________________________________________
			for ( mz=1; mz <= ms; mz++ )
			{
				//mz_ = ms - (mz - 1);ISPRAV
				N_d = N_r[mz];
				mz_ = m_t[N_d] - (mz - m_b[N_d]);
			//	if (N_d != N_r[mz_])
			//		printf ("\n OSHIBKA 5\n");
								
				for ( i=1; i <= 3; i++ )
				{
					mult_(I_s[mz],I_s_i[mz], fi3[1][i][N_d],fi3[2][i][N_d],&u,&v);
					mult_(I_s[mz_],I_s_i[mz_], Q[1][i][N_d],Q[2][i][N_d],&u1,&v1);
					A_[1][i][mz] = u + u1;
					A_[2][i][mz] = v + v1;

					mult_(I_s[mz],I_s_i[mz], Q[1][i][N_d],Q[2][i][N_d],&u,&v);
					mult_(I_s[mz_],I_s_i[mz_], fi1[1][i][N_d],fi1[2][i][N_d],&u1,&v1);
					B_[1][i][mz] = u + u1;
					B_[2][i][mz] = v + v1;
				}//i	
			}//mz

				
			// START of fill in g[mz][mz'][1:9] that does not depend on X,Y
			// see "Matrichnie elementi" p.1
			for ( mz=1; mz <= ms; mz++ )//наблюдения observation
			{
				for ( mz1=1; mz1 <= ms; mz1++ )
					for ( i=1; i <= 9; i++ )
						g_r[i][mz1] = g_i[i][mz1] = 0.0;
				
				
				//mz_ = ms - (mz - 1);ISPRAV
				N_d = N_r[mz];
				mz_ = m_t[N_d] - (mz - m_b[N_d]);
			//	if (N_d != N_r[mz_])
			//		printf ("\n OSHIBKA 1\n"); 

				k2_2 =  k2_r[N_d];//!!!!?????????
				k2_2_i = k2_i[N_d];

				f[1][1] = -alf2 + k2_2;
				f_i[1][1] = k2_2_i;
				f[2][2] = -bet2 + k2_2;
				f_i[2][2] = k2_2_i;//!!!!?????????

				gam2 = Ro - k2_2;
				gam2_i = -k2_2_i;
			
				
				for ( mz1=1; mz1 <= ms; mz1++ )//исток  source
				{
					N_d1 = N_r[mz1];
					if (N_d != N_d1)
						continue;

					if (p && q && rg_[N_d1] * fabs(Zm[mz] - Zm[mz1]) > ft_e)
						continue;
					
					//mz1_ = ms - (mz1 - 1);
					mz1_ = m_t[N_d1] - (mz1 - m_b[N_d1]);
		//			if (N_d != N_r[mz1_])
		//				printf ("\n OSHIBKA 6\n");

															
					//Неоднородная часть ФГ The inhomogeneous part of the FY					
					for ( i=0; i <= 2; i++ )
						U_c[1][i] = U_c[2][i] = 0.0;
					
					//U[][0], U[][1] & U[][2]************************************************************	

					//U[][0]=U[][01] + U[][02]; U[][2]=k^2 * U[][01] + R0 * U[][02]
														
					if (mz1 == m_t[N_d])
					{
						u = hm1[mz1] * (Wu[1][mz][mz1-1] - Wu[1][mz][mz1]);
						v = hm1[mz1] * (Wu[2][mz][mz1-1] - Wu[2][mz][mz1]);

						u1 = hm1[mz1] * (Wv[1][mz][mz1-1] - Wv[1][mz][mz1]);
						v1 = hm1[mz1] * (Wv[2][mz][mz1-1] - Wv[2][mz][mz1]);
						
						u -= I_s[mz];
						v -= I_s_i[mz];
						
						u1 -= I_c[mz];
						v1 -= I_c_i[mz];
				
					}
					else
						if (mz1 == m_b[N_d])
						{
							u = hm1[mz1+1] * (Wu[1][mz][mz1+1] - Wu[1][mz][mz1]);
							v = hm1[mz1+1] * (Wu[2][mz][mz1+1] - Wu[2][mz][mz1]);

							u1 = hm1[mz1+1] * (Wv[1][mz][mz1+1] - Wv[1][mz][mz1]);
							v1 = hm1[mz1+1] * (Wv[2][mz][mz1+1] - Wv[2][mz][mz1]);
							
							u -= I_s[mz_];
							v -= I_s_i[mz_];
							
							u1 += I_c[mz_];
							v1 += I_c_i[mz_];
						}
						else
						{
							u = hm1[mz1] * (Wu[1][mz][mz1-1] - Wu[1][mz][mz1]) +
								hm1[mz1+1] * (Wu[1][mz][mz1+1] - Wu[1][mz][mz1]);
							v = hm1[mz1] * (Wu[2][mz][mz1-1] - Wu[2][mz][mz1]) +
								hm1[mz1+1] * (Wu[2][mz][mz1+1] - Wu[2][mz][mz1]);

							u1 = hm1[mz1] * (Wv[1][mz][mz1-1] - Wv[1][mz][mz1]) +
								 hm1[mz1+1] * (Wv[1][mz][mz1+1] - Wv[1][mz][mz1]);
							v1 = hm1[mz1] * (Wv[2][mz][mz1-1] - Wv[2][mz][mz1]) +
								 hm1[mz1+1] * (Wv[2][mz][mz1+1] - Wv[2][mz][mz1]);
						}
				
					mult_(u,v, Gm1[1][N_d],Gm1[2][N_d], &u,&v);//u = U[][02]
					
					/*um = Gm1[1][N_d];
					vm = Gm1[2][N_d];
					wm = u;

					u = u * um - v * vm;
					v = v * um + wm * vm;
					*/
					U_c[1][2] = Ro * u;//R0 * U[][02]
					U_c[2][2] = Ro * v;
										
					u2 = 0.0;
					if (mz == mz1)
					{
						if (mz == m_b[N_d])
						{
							u2 = hm[mz+1] / 3.0;
							u1 -= 0.5;
						}
						else
							if (mz == m_t[N_d])
							{
								u2 = hm[mz] / 3.0;
								u1 += 0.5;
							}
							else
								u2 = (hm[mz] + hm[mz+1]) / 3.0;
							
						
					}//mz = mz1
					else
						if (mz == mz1 - 1)
						{
							u2 = hm[mz1] / 6.0;
							u1 += 0.5;
						}
						else
							if (mz == mz1 + 1)
							{
								u2 = hm[mz] / 6.0;
								u1 -= 0.5;
							}
						
					

					mult_(u1,v1, Gm1[1][N_d],Gm1[2][N_d], &u1,&v1);
					U_c[1][1] = u1;
					U_c[2][1] = v1;
					
					v2 = Gm1[2][N_d] * u2;
					u2 *= Gm1[1][N_d];//u2 = U[][01]
					
					U_c[1][0] = u + u2;
					U_c[2][0] = v + v2;	
				
					mult_(u2,v2, k2_2,k2_2_i, &u,&v);
					U_c[1][2] += u;
					U_c[2][2] += v;
						
					
					//Неоднородная часть ФГ The inhomogeneous part of the FYU[][0], U[][1] & U[][2]***********
						
				
					for ( i=1; i <= 2; i++ )//неоднородная часть The inhomogeneous part
						for ( i1=1; i1 <= 2; i1++ )
						{
							n = i + 3 * (i1 - 1);//1, 2, 4, 5
							if (i == i1)
							{
								mult_(f[i1][i],f_i[i1][i], U_c[1][0],U_c[2][0], &u,&v);
								g_r[n][mz1] = u;
								g_i[n][mz1] = v;
							}
							else
							{
								g_r[n][mz1] = f[i1][i] * U_c[1][0];
								g_i[n][mz1] = f[i1][i] * U_c[2][0];						
							}
						}

						//???????????????????????!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!				
					
					
					g_r[9][mz1] =  U_c[1][2];
					g_i[9][mz1] =  U_c[2][2];

											
					for ( i1=0, n=6, i=1; i <= 2; i++ )
					{
						i1 +=3;//3, 6
						n++;//7, 8 
						g_r[i1][mz1] = f[i][3] * U_c[1][1];//!!//
						
						g_r[n][mz1] = f[3][i] * U_c[1][1];
						
						g_i[i1][mz1] = f[i][3] * U_c[2][1];//!!!//
						g_i[n][mz1] = f[3][i] * U_c[2][1];	
						
					}//n

					
				}//mz1 - неоднородная часть The inhomogeneous part END
				
								
				//if (buv)// odnorodnaia 
				//однородная часть ФГ homogeneous part of the FY
				SR = 1.0;
				for ( mz1=1; mz1 <= ms; mz1++ )//исток  source
				{
					
					if (buv)// odnorodnaia 
					{	

						//mz1_ = ms - (mz1 - 1);
						N_d1 = N_r[mz1];
						mz1_ = m_t[N_d1] - (mz1 - m_b[N_d1]);
			//			if (N_d1 != N_r[mz1_])
			//				printf ("\n OSHIBKA 9\n");
					
						//___________________________________//  V1,2,3,4,5  without  delta
						//bur2 = bur1 = (ro_ * (Zm[mz] - Zm[mz1]) < Rex) ? 1 : 0;//?????????????????
						
						bur2 = bur1 = 1;
					//	if (!bur1)
					//		continue;
//___________________________________________________________________________________________________________________
						if (N_d == N_d1)
						{
							for ( i=1; i <= 3; i++ )
							{							
								u = v = u1 = v1 = 0.0;
								//if (bur1)
									mult_(I_s[mz],I_s_i[mz], A_[1][i][mz1],A_[2][i][mz1], &u,&v);
								//if (bur2)
									mult_(I_s[mz_],I_s_i[mz_], B_[1][i][mz1],B_[2][i][mz1], &u1,&v1);
														
								VE[1][i] = u + u1;
								VE[2][i] = v + v1;
							
								if ( i > 1 )
								{								
									u = v = u1 = v1 = 0.0;
								//	if (bur1)
										mult_(I_c[mz],I_c_i[mz], A_[1][i][mz1],A_[2][i][mz1], &u,&v);
								//	if (bur2)
										mult_(I_c[mz_],I_c_i[mz_], B_[1][i][mz1],B_[2][i][mz1], &u1,&v1);
																				
									VE[1][i+2] = u - u1;
									VE[2][i+2] = v - v1;						
								}
							}//i

							
						}//N_d == N_d1
						else
						//	if (abs(N_d - N_d1) > 10)
						//		continue;
						//	else
						{
							divide_(EPS_r[N_d1],EPS_i[N_d1], EPS_r[N_d],EPS_i[N_d], &C_[0],&D_[0]);
							for ( i=1; i <= 3; i++ )
							{									
								if (N_d < N_d1 )//(.) наблюдения выше (.) истока M1
								{
									D_1[1] = D_t[1][i][N_d-1][N_d1];
									D_1[2] = D_t[2][i][N_d-1][N_d1];

									D_[1] = D_t[1][i][N_d][N_d1];
									D_[2] = D_t[2][i][N_d][N_d1];

									C_[1] = A_[1][i][mz1];
									C_[2] = A_[2][i][mz1];
								}
								else//(.) наблюдения ниже (.) истока M3
								{
									D_1[1] = D_b[1][i][N_d-1][N_d1];
									D_1[2] = D_b[2][i][N_d-1][N_d1];

									D_[1] = D_b[1][i][N_d][N_d1];
									D_[2] = D_b[2][i][N_d][N_d1];

									C_[1] = B_[1][i][mz1];
									C_[2] = B_[2][i][mz1];
								}
								
								//if (bur1)
								{
									mult_(D_1[1],D_1[2], I_s[mz],I_s_i[mz],&M_1,&M_11);
									mult_(D_[1],D_[2], I_s[mz_],I_s_i[mz_],&u,&v);
									mult_(M_1 + u,M_11 + v, C_[1],C_[2], &u,&v);

									if (i == 3)//
										mult_(u,v, C_[0],D_[0], &u,&v);
							
								}
								//else
								//	u = v = 0.0;
														
								VE[1][i] = u;
								VE[2][i] = v;
							
								if ( i > 1 )
								{								
									//if (bur1)
									{
										mult_(D_1[1],D_1[2], I_c[mz],I_c_i[mz],&M_1,&M_11);
										mult_(D_[1],D_[2], I_c[mz_],I_c_i[mz_],&u,&v);
										mult_(M_1 - u,M_11 - v,C_[1],C_[2], &u,&v);

										if (i == 3)
											mult_(u,v, C_[0],D_[0], &u,&v);
										else
										{
											mult_(u,v, G2[1][N_d1],G2[2][N_d1], &u,&v);
											divide_(u,v, G2[1][N_d],G2[2][N_d], &u,&v);
										}
									}
									//else
									//	u = v = 0.0;
																		
									VE[1][i+2] = u;
									VE[2][i+2] = v;						
								}//i > 1
																	
							}//i

												
						}//N_d != N_d1
						
						
//____________________________________________________________________________________________________________________
						{
							double u11,v11,v12,u12;
							mult_(VE[1][1],VE[2][1], k2_r[N_d1],k2_i[N_d1], &u11,&v11);
						
							mult_(VE[1][2],VE[2][2], G2[1][N_d1],G2[2][N_d1], &u12,&v12);
							
							//g_r[1][mz][mz1] += br * u11 - ar * u12;
							//g_i[1][mz][mz1] += br * v11 - ar * v12;
							G_r[1] = g_r[1][mz1] + br * u11 - ar * u12;
							G_i[1] = g_i[1][mz1] + br * v11 - ar * v12;

							//if (abs(p) > p_g  && abs(q) > q_g && fabs(G_r[1] / SR) <= ft_l)
							if (p  && q && fabs(G_r[1] / SR) <= ft_l)
						//	if (p && q &&  fabs(rg_[N_d] * Zm[mz] - rg_[N_d1] *Zm[mz1]) > ft_e)
								continue;

							u = -abr * (u11 + u12);
							G_r[2] = g_r[2][mz1] + u;//
							G_r[4] = g_r[4][mz1] + u;
							//g_r[2][mz][mz1] += u;//
							//g_r[4][mz][mz1] += u;
					
							u = -abr * (v11 + v12);
							//g_i[2][mz][mz1] += u;//
							//g_i[4][mz][mz1] += u;
							G_i[2] = g_i[2][mz1] + u;//
							G_i[4] = g_i[4][mz1] + u;
					
					
							//g_r[5][mz][mz1] += ar * u11 - br * u12;
							//g_i[5][mz][mz1] += ar * v11 - br * v12;
							G_r[5] = g_r[5][mz1] + ar * u11 - br * u12;
							G_i[5] = g_i[5][mz1] + ar * v11 - br * v12;
										
							//g_r[9][mz][mz1] += Ro * VE[1][3];
							//g_i[9][mz][mz1] += Ro * VE[2][3];
							G_r[9] = g_r[9][mz1] + Ro * VE[1][3];
							G_i[9] = g_i[9][mz1] + Ro * VE[2][3];

							for ( i1=0, n=6, i=1; i <= 2; i++ )
							{
								i1 +=3;//3, 6
								n++;//7, 8 
								//g_r[i1][mz][mz1] += f[i][3] * VE[1][5];
								//g_i[i1][mz][mz1] += f[i][3] * VE[2][5];
								G_r[i1] = g_r[i1][mz1] + f[i][3] * VE[1][5];
								G_i[i1] = g_i[i1][mz1] + f[i][3] * VE[2][5];
						
								//g_r[n][mz][mz1] += f[3][i] * VE[1][4];
								//g_i[n][mz][mz1] += f[3][i] * VE[2][4];
								G_r[n] = g_r[n][mz1] + f[3][i] * VE[1][4];
								G_i[n] = g_i[n][mz1] + f[3][i] * VE[2][4];
						
							}
						}
					
						for ( i=1; i <= 9; i++ )
						{
							G_r[i] *= 0.5;
							G_i[i] *= 0.5;
						}
					}// buv
					else
					{
						for ( i=1; i <= 9; i++ )
						{
							G_r[i] = g_r[i][mz1] * 0.5;
							G_i[i] = g_i[i][mz1] * 0.5;
						}
					}//!buv	


					//if (fabs(G_r[1] / SR) <= ft_l)
					//if (abs(p) > p_g  && abs(q) > q_g && fabs(G_r[1] / SR) <= ft_l)
					if (p  && q && fabs(G_r[1] / SR) <= ft_l)
				//	if (p && q &&  fabs(rg_[N_d] * Zm[mz] - rg_[N_d1] * Zm[mz1]) > ft_e)
						continue;
								
					MATR_(mz,mz1);

					//if (!p && !q && mz == mz1 == 1)
					//if (!p && !q && mz == 1 && mz1 == 1)
					//if (!p && !q && mz + mz1 == 2)
					if (p && q && mz + mz1 == 2)
						SR = g_r[1][1];

					}//mz1 /odn
						
			}//mz
			
			
			
		}//q
	}//p

	
	if (Ekv)
		DIAGON(); 
	else
		DIAGON_N(); 

	
	CONVERG( );//улучшение сходимости improvement of convergence//	// asymptotic calc. of integrals outside diffr.scheme


	if (N_ell > 1)
	{
		for ( m=0; m < MS_; m++ )
			for ( m1=0; m1 < MS_; m1++ )
			{
				u = Alap[m][m1].real;
				v = Alap[m][m1].imag;

				for ( i=2; i <= N_ell; i++ )
				{
					m_ = m + MS_* (i - 1);
					m1_ = m1 + MS_* (i - 1);

					Alap[m_][m1_].real = u;
					Alap[m_][m1_].imag = v;

				}//i
			}//m1
	}//N_ell > 1
	

	
	CLineEqC(Alap,MS_N,Blap,N_corn);
	// solves lin.eq.system using MKL; result in Blap

	
	for ( Pol=0; Pol < 2; Pol++ )
	{
	//	u = v = 0.0;
	for ( m=0; m < MS_N; m++ )		
		{
			ASR_[Pol][m+1] = Blap[Pol][m].real;
			ASI_[Pol][m+1] = Blap[Pol][m].imag;	
	//		u += ASR_[Pol][m+1];
	//		v += ASI_[Pol][m+1];
			
		}
	//	printf ("\n	     new 2   	  %7.5e       %7.5e", u,v );
	//	fprintf( out,"\n		    new  2   	  %7.5e       %7.5e", u,v );
	}

	// X(i,mr,mf,mz)  [see "Presentation" eq(2)] is the solution of lin.eq.systemis in ASR_ and ASI_
	// i=1,2,3 for Ex,Ey,Ez 

	D_N_PER( ); //рассчитывается поле в дальней зоне решетки calculated field in the far zone of the lattice

	FM6();
	
	fmat_2d( (void**)rte,1,3, 1, sizeof( double) );
	fmat_2d( (void**)ite,1,3, 1, sizeof( double) );
	fmat_2d( (void**)rq,1,3, 1, sizeof( double) );
	fmat_2d( (void**)iq,1,3, 1, sizeof( double) );
	fmat_2d( (void**)rq_,1,3, 1, sizeof( double) );
	fmat_2d( (void**)iq_,1,3, 1, sizeof( double) );

	fmat_3d( (void***)A_, 1,2,1,3, 1, sizeof( double) );
	fmat_3d( (void***)B_, 1,2,1,3, 1, sizeof( double) );

	
}/* F */

