
#include "VIBRr.h"

//double Rek, Imk, ky, r_tau[3], i_tau[3], r_mukap, i_mukap;


void F_df_IE_MK(int PP, double *R, double *T)//конечный
{
	double  kk, **RX, **IX, x, y, u, v,  r_Matr, i_Matr,t, rR, iR, rT, iT,rU,iU,udd,vdd;
	double ky,**RA, **IA, alf, gam, gam2,    u_m_mu,c,s,Jm,al_L,rY,iY,zz,rUm=0.0,iUm=0.0,
		kap_1[4], kap[4], xa, ya, rV, iV, rro, iro, rTA, iTA,ru,iu,rex,iex,
		rR_ext, iR_ext, rT_ext, iT_ext,cn,cnu;
	int i, j, m, n, m_u, nu, p, cmu;
	int vardd = 1, N_g_;
	
	N_g_ = 3 * N_g;

	if (kbhit())
		if (getch() == ESCape)
		{
			fflush(out);
			exit(1);
		}
	
	L = 0.5 * l;
	H = 0.5 * h;

	RA = (double**)amat_2d(1, MS, 1, MS + 1, sizeof(double));
	IA = (double**)amat_2d(1, MS, 1, MS + 1, sizeof(double));
	RX = (double**)amat_2d(1, M_X, 1, M_Y, sizeof(double));
	IX = (double**)amat_2d(1, M_X, 1, M_Y, sizeof(double));
	
	ky = k * cos(TET);
	t = 0.5 / d;
	
//________________________________________
	kk = kx * kx;
	
	/*mult(EPS[1], EPS[2], MU[1], MU[2], &x, &y);

	r_tau[3] = k2 *  (x - 1.0);
	i_tau[3] = k2 * y;
	*/
	epsG[2] = eps[2];
	epsG[3] = eps[3];
	muG[2] = mu__[2];
	muG[3] = mu__[3];
	h2G = h2;

	E_ext(PP, 0, &rR_ext, &iR_ext, &rT_ext, &iT_ext);

	if (TEST == 3)
	{
		E_ext(PP, 1, &rR_ext, &iR_ext, &rT_ext, &iT_ext);//!!!!!!!!!!!!!!!
		u = rR_ext * rR_ext + iR_ext * iR_ext;
		v = rT_ext * rT_ext + iT_ext * iT_ext;

		printf("\n s(anal)  f=%6.4lf   R=%6.4e    T=%6.4e    P=%6.4e ", f, u, v, 1.0 - u - v);
		//exit(1);
	}

	if (PP)
	{
		r_tau[1] = k2 *  (MU[1] - 1.0);
		i_tau[1] = k2 * MU[2];

		divide(EPS[1] - 1.0, EPS[2], EPS[1], EPS[2], &u, &v);
		r_tau[2] = u;
		i_tau[2] = v;

		divide(1.0, 0.0, EPS[1], EPS[2], &u, &v);
		r_tau[3] = k2 *  (MU[1] - u);
		i_tau[3] = k2 * (MU[2] - v);

		divide(1.0, 0.0, EPS[1], EPS[2], &u, &v);

		if (vardd)
		{
			udd = u;
			vdd = v;
			zz = 1.0;
		}
		else
		{
			udd = 1.0;
			vdd = 0.0;
			zz = 1.0;
		}
	}
	else
	{
		r_tau[1] = k2 *  (EPS[1] - 1.0);
		i_tau[1] = k2 * EPS[2];

		divide(MU[1] - 1.0, MU[2], MU[1], MU[2], &u, &v);
		r_tau[2] = u;
		i_tau[2] = v;
		
		divide(1.0, 0.0, MU[1], MU[2],  &u, &v);
		if (vardd)
		{
			udd = u;
			vdd = v;
						
			zz = 1.0;
		}
		else
		{
			udd = 1.0;
			vdd = 0.0;
			zz = 1.0;
		}
		r_tau[3] = k2 *  (EPS[1] - u);
		i_tau[3] = k2 * (EPS[2] - v);
		
	}
	
	RA = (double**)amat_2d(1, MS, 1, MS + 1, sizeof(double));
	IA = (double**)amat_2d(1, MS, 1, MS + 1, sizeof(double));
	
	for (m = 1; m <= MS + 1; m++)
	for (n = 1; n <= MS; n++)
		RA[n][m] = IA[n][m]= 0.0;

	for (j = 0, m = 1; m <= M_X; m++)
	for (cn = -1.0, n = 1; n <= M_Y; n++)
	{
		j++;

		cn = -cn;
		
		u = 1.0 / (L * Pi);
		u /= (1.0 + 2.0 *(double)(m - 1));
		v = 1.0 / (H * Pi);
		v /= (1.0 + 2.0 *(double)(n - 1));
		u *= v;
		
		
		RA[j][j] = u * udd;
		IA[j][j] = u * vdd;
		
		y = ky * H;

		u_m_mu = J(m - 1, kx * L) * J(n - 1, y);

		x = cos(y);
		y = sin(y);
		mult(x, -y, rR_ext, iR_ext, &u, &v);

	//	rc = cn;
		RA[j][MS + 1] = u_m_mu * (x + u * cn);
		IA[j][MS + 1] = u_m_mu * (y + v * cn);
		
		for (i = 0, cmu = -1, m_u = 1; m_u <= M_X; m_u++)
		{
			cmu = -cmu;
			for (cnu = -1.0, nu = 1; nu <= M_Y; nu++)
			{
				i++;
				cnu = -cnu;
								
				r_Matr = i_Matr = 0.0;
				for (p = -N_g_; p <= N_g_; p++)
				{
					alf = Pi_2 / d * (double)p + kx;

					F_GRIN_dop(PP, alf, &rro, &iro, &rTA, &iTA);

					al_L = alf * L;

					gam2 = alf * alf - k2;
					
					// _____________________________________A .	поверхностные интегралы
					Jm = J(m - 1, al_L);
					u_m_mu = Jm *  J(m_u - 1, al_L);
					
					
					U_(n - 1, nu - 1, (int)cn, (int)cnu, gam2, &rU, &iU);
					mult(r_tau[3], i_tau[3], rU, iU, &ru, &iu);// особая

					V_(PP, n, nu, alf, &rV, &iV);
					
					u = alf * alf + gam2 /** (cn)*/;
					r_tau[4] = r_tau[1] + u * r_tau[2];
					i_tau[4] = i_tau[1] + u * i_tau[2];
					mult(r_tau[4], i_tau[4], rV, iV, &x, &y);
					mult(x, y, rro, iro, &x, &y);

					xa = u_m_mu * (ru + x);
					ya = u_m_mu * (iu + y);

					r_Matr += xa;
					i_Matr += ya;
										
					
					// _____________________________________A end a.	поверхностные интегралы
					s = -alf * Jm * Q_mu(m_u - 1, L, -1, al_L);
					
					mult(rV, iV, rro, iro, &x, &y);

					u = s * (rU + x);
					v = s * (iU + y);
					mult(u, v, r_tau[2], i_tau[2], &x, &y);

					r_Matr += x;
					i_Matr += y;
					// _________________________________________b.	END Контурные интегралы B

					//if (cn == cnu)
					{
							s = -u_m_mu / (H * sqrt(Pi_2));

						//Y
						if (gam2 >= 0.0)
						{
							gam = sqrt(gam2);
							x = I_ex(n-1,gam * H);
							
							//c = 0.5 * Pi * (double)(nu - n);
							c = 0.5 * Pi * (double)(n - 1);
							iY = sin(c) * x;
							rY = cos(c) * x;

							rex = (exp(-2.0 * gam * H) - cnu)/* * cn*/;
							iex = 0.0;

							c = 0.5 * Pi * (double)(nu - 1);
							/*si = sin(c);
							co = cos(c);*/
							
						}//gam2 >= 0.0
						else
						{
							gam = sqrt(-gam2);
																					
							x = J(n - 1, gam * H);

							y = -gam * H + Pi * (double)(n - 1);
							iY = sin(y) * x;
							rY = cos(y) * x;

							y = -2.0 * gam * H;
							rex = (cos(y) - cnu)/* * cn*/;
							iex = sin(y)/* * cn*/;
						}//gam2 < 0.0
						//Y END
						
						mult(-rro, -iro, rex, iex, &rex, &iex);
					
						rex += cn + cnu;
						
						mult(s * rY, s * iY, r_tau[2], i_tau[2], &x, &y);

						mult(x, y, rex, iex, &x, &y);

						v = Pi05 * (double)(nu - 1);
						u = cos(v);
						v = sin(v);
						mult(x, y, u,v, &x, &y);

						r_Matr += x;
						i_Matr += y;
					}//cn != cnu
					// _________________________________________END Контурные интегралы W2

				/*	F_df_dop( PP,  m,  m_u,  n,  nu,  alf,  &x,&y);
					r_Matr += x;
					i_Matr += y;*/
				}//p

				RA[j][i] -= r_Matr * t;
				IA[j][i] -= i_Matr * t;

			
			}//nu
			
		}//mu
	}//j

	

	CELINEG(1, MS, RA, IA);
	//________________________________________//Решение СЛАУ

	
	
	for (j=0, m = 1; m <= M_X; m++)
	for (n = 1; n <= M_Y; n++)
	{
		j++;
		RX[m][n] = RA[j][MS+1];
		IX[m][n] = IA[j][MS+1];
		//printf("\n      %2d  %2d         RX=%6.4e    IX=%6.4e", n, m, fabs(RX[m][n]), fabs(IX[m][n]));
	}
	

	fmat_2d((void**)RA, 1, MS, 1, sizeof(double));
	fmat_2d((void**)IA, 1, MS, 1, sizeof(double));

	//________________________________________//Дифрагированное поле
	{
		double rr, ir, rt, it, t, si, co, jal, jgam,  gH,   kap_, u_sc, v_sc, u_sc2, v_sc2,
			rYp, iYp, rYm, iYm;

		r_Matr = i_Matr = 0.0;
		for (p = -N_g; p <= N_g; p++)
		{
			alf = Pi_2 / d * (double)p + kx;

			gam2 = alf * alf - k2;
			if (gam2 >= 0.0)
			{
				if (p < 0)
					continue;
				else
					break;
			}
			
			gam = sqrt(-gam2);
			t = 0.5 / (d * gam);
			
			gH = gam * H;
		
			co = cos(gH);
			si = sin(gH);
						
			kap_ = sqrt(-alf * alf + k2 * eps[3] * mu__[3]);
			
			F_GRIN_dop(PP, alf, &rro, &iro, &rTA, &iTA);

			if (!p)
			{
				rR = rR_ext;
				iR = iR_ext;
				rT = rT_ext;
				iT = iT_ext;

				/*u = rR * rR + rT * rT + iR * iR + iT * iT;
				printf("\n  PPPPPPPPPP=%6.4e", 1.0 - u);*/
			}
			else
				rR = iR = rT = iT = 0.0;

						
			for (cmu = -1, m_u = 1; m_u <= M_X; m_u++)
			{	
				cmu = -cmu;

				jal = J(m_u - 1,  alf * L);
			
										
				for (cnu = -1, nu = 1; nu <= M_Y; nu++)
				{		
					cnu = -cnu;
												
					//			R1 & T1______________________________________________	
					jgam = J(nu - 1, gH);
					
					u = jal *  jgam;
					x = u * r_tau[1];
					y = u * i_tau[1];

					mult(co, -si, rro, iro, &u_sc, &v_sc);
					u_sc += (double)cnu * co;
					v_sc += (double)cnu * si;
					mult(x, y, u_sc, v_sc, &rr, &ir);

					mult(co, -si, rTA, iTA, &u_sc2, &v_sc2);
					mult(x, y, u_sc2, v_sc2, &rt, &it);
				
					//			R1 & T1  END______________________________________________	
				
					x = -alf * jgam * (Q_mu(m_u - 1, L, -1, alf  * L) - alf * jal);
					y = i_tau[2] * x;
					x *= r_tau[2];

					mult(x, y, u_sc, v_sc, &u, &v);
					rr += u;
					ir += v;
					
					mult(x, y, u_sc2, v_sc2, &u, &v);
					rt += u;
					it += v;

					//			R2 & T2  END______________________________________________	

					YSH(nu, gam, &rYp, &iYp);
					YSH(nu, -gam, &rYm, &iYm);
					/*rYm = cnu * rYp;
					iYm = cnu * iYp;*/
					
					x = jal * gam;
					y = i_tau[2] * x;
					x *= r_tau[2];

					mult(rYp, iYp, rro, iro, &u, &v);
					u += rYm;
					v += iYm;
					mult(x, y, u, v, &u, &v);
					rr += u;
					ir += v;
					

					x = -jal * kap_;
					y = i_tau[2] * x;
					x *= r_tau[2];
					mult(rYp, iYp, rTA, iTA, &u, &v);
					mult(x, y, u, v, &u, &v);
					rt += u;
					it += v;

					//			R3 & T3  END______________________________________________	
					
					//__________________________________________________________________________________

					mult(RX[m_u][nu], IX[m_u][nu], rr, ir, &rr, &ir);
					mult(RX[m_u][nu], IX[m_u][nu], rt, it, &rt, &it);
									
					rR += t * ir;
					iR -= t * rr;

					rT += t * it;
					iT -= t * rt;

				}//nu
			}//mu

			
			rR = rR * rR + iR * iR;

			kap[3] = sqrt(k2 * eps[3] * mu__[3] - kx * kx);

			x = (PP) ? eps[3] : mu__[3];
			kap_1[3] = kap[3] / x;
			rT = (rT * rT + iT * iT) * kap_1[3] / gam;
			printf("\n p=%2d      R=%6.4e    T=%6.4e", p, rR,rT);

			r_Matr += rR;
			i_Matr += rT;
		}//p

		*R = r_Matr;
		*T = i_Matr;

	}	//________________________________________//Дифрагированное поле
	
}//F_df

//u = 0.0;
//for (m = 1; m <= MS; m++)
//{
//	fprintf(out, "\n");
//	for (n = 1; n <= m; n++)
//	{
//		u += RA[n][m]+ IA[n][m];
//		if (n == m)
//		printf("\n  %2d  %2d         R=%6.4e    T=%6.4e", n, m, RA[n][m], IA[n][m]);
//		fprintf(out, "\n  %2d  %2d       %6.4e    %6.4e        %6.4e    %6.4e",
//			n, m, RA[n][m], IA[n][m], RA[m][n], IA[m][n]);
//		
//	}
//}
//printf("\n  u=%6.4e", u);
//printf("\n  u=%6.4e  u=%6.4e", rUm,iUm);