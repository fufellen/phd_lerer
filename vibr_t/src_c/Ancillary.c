
#include "VIBRr.h"


void Ancillary( void)// 
{ // !!!core function 
	int N_d,mz,mz1,mz_,mz1_,bu1,buz;
	double u,u1,v,v1,u2,v2,x,y,sh_[5][221],sz[3],cz[3],sp[5],cp[5],cm[5],sm[5],Z,Z1,Z_,Z1_;
	
	for ( N_d=2; N_d <= N_hol + 1; N_d++ )
	{
		x = rg_[N_d] * H_[N_d-1];
		y = ig_[N_d] * H_[N_d-1];
		C_sh(x,y, &u,&v);
		divide_(1.0,0.0, u,v, &sh_[1][N_d],&sh_[2][N_d]);//1.0 /sinh(gamma*h)

		mult_(rg_[N_d],ig_[N_d], u,v, &u,&v);
		divide_(1.0,0.0, u,v, &sh_[3][N_d],&sh_[4][N_d]);//1.0 / (gamma * sinh(gamma*h))
	}
			
	for ( mz=1; mz <= ms; mz++ )//calculate sinh(gamma*Zm),  cosh(gamma*Zm)
	{ // see "Matrichnie elementi" p4
		N_d = N_r[mz];
		Cn_sh_ch(rg_[N_d],ig_[N_d],Zm[mz] - h_b[N_d],H_[N_d-1], &u1,&v1, &u2,&v2);

		mult_(u1,v1, sh_[1][N_d] ,sh_[2][N_d] , &u1,&v1);
		Sh[1][mz] = u1;//sinh(gamma*Zm) / sinh(gamma*H)
		Sh[2][mz] = v1;
		
		mult_(u2,v2, sh_[1][N_d] ,sh_[2][N_d] , &u1,&v1);
		Ch[1][mz] = u1;//cosh(gamma*Zm) / sinh(gamma*H)
		Ch[2][mz] = v1;
				
	}//mz

	for ( mz=1; mz <= ms; mz++ )//calculateÂà÷èñëåíèåI_s[mz],I_c[mz]
	{ 
		N_d = N_r[mz];
		INT(N_d,1,mz, rg_[N_d], ig_[N_d], Sh,Ch, sz,cz);
		
		I_s[mz] = sz[1];
		I_s_i[mz] = sz[2];

		mult_(rg_[N_d], ig_[N_d],cz[1],cz[2], &u,&v);
		I_c[mz] = u;
		I_c_i[mz] = v;
		
	}//mz
	
	for ( mz=1; mz <= ms; mz++ )//calculate U(Zm,Zk),  V(Zm,Zk), V(Zm-0,Zm)
	{
		//mz_ = ms - (mz - 1);ISPRAV
		
		N_d = N_r[mz];
		mz_ = m_t[N_d] - (mz - m_b[N_d]);
	//	if (N_d != N_r[mz_])
	//				printf ("\n OSHIBKA 2\n");

		Z = Zm[mz] - h_b[N_d];
		Z_ = Zm[mz_] - h_b[N_d];
			
		for ( mz1=1; mz1 <= ms; mz1++ )//mz >= mz1
		{
			if (N_d != N_r[mz1])
				continue;
			
			//mz1_ = ms - (mz1 - 1);ISPRAV
			mz1_ = m_t[N_d] - (mz1 - m_b[N_d]);

			Z1 = Zm[mz1] - h_b[N_d];
			Z1_ = Zm[mz1_] - h_b[N_d];
		

		//	if (N_d != N_r[mz1_])
			//		printf ("\n OSHIBKA 3\n");
			
			Cn_sh_ch(rg_[N_d],ig_[N_d],Z_ - Z1,H_[N_d-1], &sm[1],&sm[2], &cm[1],&cm[2]);

			if ( mz >= mz1)
			{
				Cn_sh_ch(rg_[N_d],ig_[N_d],Z_ + Z1,H_[N_d-1], &sp[1],&sp[2], &cp[1],&cp[2]);
				
				sp[1] = -sp[1];//+0
				sp[2] = -sp[2];

				if (mz == mz1)
				{
					sp[3] = -sp[1];//-0
					sp[4] = -sp[2];
				}
			}
			else
				Cn_sh_ch(rg_[N_d],ig_[N_d],Z + Z1_,H_[N_d-1], &sp[1],&sp[2], &cp[1],&cp[2]);
				
			u1 = cp[1] - cm[1];
			v1 = cp[2] - cm[2];
			mult_(u1,v1, sh_[3][N_d], sh_[4][N_d], &u1,&v1);
			U[1][mz][mz1] = u1 * 0.5;
			U[2][mz][mz1] = v1 * 0.5;
			
			u2 = sp[1] + sm[1];
			v2 = sp[2] + sm[2];
			mult_(u2,v2, sh_[1][N_d],sh_[2][N_d], &u1,&v1);
			V[1][mz][mz1] = u1 * 0.5;
			V[2][mz][mz1] = v1 * 0.5;
					
			if (mz == mz1)
			{
				u2 = sp[3] + sm[1];
				v2 = sp[4] + sm[2];
			
				mult_(u2,v2, sh_[1][N_d],sh_[2][N_d], &u1,&v1);
							
				Vd[1][mz] = u1 * 0.5;
				Vd[2][mz] = v1 * 0.5;
			}
		}//mz1
	}//mz

	
	for ( mz=1; mz <= ms; mz++ )//calculate   W(Zm,Zk)
	{
		N_d = N_r[mz];

		//bu1 = (mz == 1) ? 0 : 1; 
		//buz = (mz == ms) ? 0 : 1; 
		bu1 = (mz == m_b[N_d]) ? 0 : 1; 
		buz = (mz == m_t[N_d]) ? 0 : 1; 
			
			
		for ( mz1=1; mz1 <= ms; mz1++ )//mz >= mz1
		{
			//mz1_ = ms - (mz1 - 1);ISPRAV

			if (N_d != N_r[mz1])
				continue;

			mz1_ = m_t[N_d] - (mz1 - m_b[N_d]);
		//	if (N_d != N_r[mz1_])
		//			printf ("\n OSHIBKA 7\n");
			//_____________________________Wu
				//_____________________Wu(2)
				
				u = v = 0.0;
				if (bu1)
				{
					u = hm1[mz] * (U[1][mz-1][mz1] - U[1][mz][mz1]);
					v = hm1[mz] * (U[2][mz-1][mz1] - U[2][mz][mz1]);
					
				}

				if (buz)
				{
					u += hm1[mz+1] * (U[1][mz+1][mz1] - U[1][mz][mz1]);
					v += hm1[mz+1] * (U[2][mz+1][mz1] - U[2][mz][mz1]);
				}
							
				//___________________End Wu(2)
				

				//_____________________Wu(1)

				
				if (!bu1)
				{
					u -= Sh[1][mz1_];
					v -= Sh[2][mz1_];
				}
				else
					if (!buz)
					{
						u -= Sh[1][mz1];
						v -= Sh[2][mz1];
					}
					
					if (mz == mz1)
						u += 1.0;
					
			
			//_____________________End Wu(1)
			
				mult_(u,v, Gm1[1][N_d],Gm1[2][N_d], &u1,&v1);
			
				Wu[1][mz][mz1] = u1;
				Wu[2][mz][mz1] = v1;
			
			//____________________________________End Wu

			if (mz == mz1)
			{
				u1 = Vd[1][mz];
				v1 = Vd[2][mz];
			}
			else
			{
				u1 = V[1][mz][mz1];
				v1 = V[2][mz][mz1];
			}

			if (mz + 1 == mz1)
			{
				u2 = Vd[1][mz1];
				v2 = Vd[2][mz1];
			}
			else
			{
				u2 = V[1][mz+1][mz1];
				v2 = V[2][mz+1][mz1];
			}

			u = v = 0.0;
			if (bu1)
			{
				u = hm1[mz] * (V[1][mz-1][mz1] - u1);
				v = hm1[mz] * (V[2][mz-1][mz1] - v1);
			}
			
			if (buz)
			{
				u += hm1[mz+1] * (u2 - V[1][mz][mz1]);
				v += hm1[mz+1] * (v2 - V[2][mz][mz1]);
			}
			
			mult_(u,v, Gm1[1][N_d],Gm1[2][N_d], &u,&v);
			
			Wv[1][mz][mz1] = u;
			Wv[2][mz][mz1] = v;
		}//mz1
	}//mz

			
}//Ancillary


void INT(int N_d,int h,int m, double kap_r, double kap_i, double **S, double **C, double *s, double *c)
{
	double xc,yc, u,v,u1,v1,uc,vc,uc1,vc1;

	mult_(kap_r,kap_i,kap_r,kap_i, &xc,&yc);
	
	if (h)//hiperbol
	{
		xc = - xc;
		yc = - yc;
	}

	if (m != m_b[N_d])
	{
		divide_(S[1][m] - S[1][m-1],S[2][m] - S[2][m-1], xc,yc, &u,&v);
		u *= hm1[m];
		v *= hm1[m];

		divide_(C[1][m] - C[1][m-1],C[2][m] - C[2][m-1], xc,yc, &uc,&vc);
		uc *= hm1[m];
		vc *= hm1[m];
	}//m != 1
	else
	{
		divide_(C[1][m],C[2][m], kap_r, kap_i, &u,&v);
		if (h)
		{
			u = - u;
			v = - v;
		}
			
		divide_(-S[1][m],-S[2][m], kap_r, kap_i, &uc,&vc);
		
	}//m == 1

	if (m != m_t[N_d])
	{
		divide_(S[1][m] - S[1][m+1],S[2][m] - S[2][m+1], xc,yc, &u1,&v1);
		u1 *= hm1[m+1];
		v1 *= hm1[m+1];
		
		divide_(C[1][m] - C[1][m+1],C[2][m] - C[2][m+1], xc,yc, &uc1,&vc1);
		uc1 *= hm1[m+1];
		vc1 *= hm1[m+1];
	}//m != M_z
	else
	{
		divide_(-C[1][m],-C[2][m], kap_r, kap_i, &u1,&v1);
		if (h)//hiperbol
		{
			u1 = - u1;
			v1 = - v1;
		}
		
		divide_(S[1][m],S[2][m], kap_r, kap_i, &uc1,&vc1);
	}//m == M_z

	s[1] = u1 + u;
	s[2] = v1 + v;

	c[1] = uc1 + uc;
	c[2] = vc1 + vc;

}//INT