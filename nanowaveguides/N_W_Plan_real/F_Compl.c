
#include "VIBRr.h"


void F_Compl( double *KZ,double *F_C) 
{ // !!!core function 
	//  
	int m;
	double alf_[3],alf2_[3],u2,v2,u1,v1,u,v,rk,ik,rk_[111],ik_[111],us,vs, uc,vc,ch[3],sh[3],
		Q[111],T[111], Qi[111],Ti[111], D[111],Di[111];
	
	
	if( kbhit() )
	   if( getch() == ESCape )
		 {
		 fflush( out );
		 exit( 1 );
		 }
	   	
	alf_[1] = k * KZ[1];// комплексная постоянная распространения 
	alf_[2] = -k * KZ[2];
	mult_a(alf_,alf_,alf2_);

	for ( m=1; m <= Nsl + 1; m++ )
	{
		k2_r[m] = k2 * EPS_r[m];
		k2_i[m] = k2 * EPS_i[m];

		u1 = alf2_[1] - k2_r[m]; // GammaSquared in each layer//ISPR!!!!!!!!!!!!!!!!!!!!!!!!!
		v1 =  alf2_[2] - k2_i[m];
		sq_compl(u1,v1, &u,&v);

		rk = u;//kapa
		ik = v;

		if (E)
			divide(u,v,  EPS_r[m],EPS_i[m],  &u,&v);

		rk_[m] = u;//kapa/eps
		ik_[m] = v;
		
		if ( H_[m-1] > 1.0e-7)
		{
			rk *= H_[m-1];
			ik *= H_[m-1];

			if (Nsl <= 2 )
			{
				sinh_C(rk,ik, &us,&vs);
				sh[1] = us;//sinh(kapa*h)
				sh[2] = vs;

				cosh_C(rk,ik, &uc,&vc);
				ch[1] = uc;//cosh(kapa*h)
				ch[2] = vc;
			}//Nsl <= 2
			else//Nsl > 2
			{
				SH_CTH(rk,ik, &us,&vs,&uc, &vc);

				mult(rk_[m],ik_[m],us,vs,&us,&vs);
				Q[m] = us;
				Qi[m] = vs;

				mult(rk_[m],ik_[m],uc,vc,&uc,&vc);
				T[m] = uc;
				Ti[m] = vc;
			}//Nsl > 2
		}//H_[m-1] > 1.0e-7
		else
		{
			T[m] = rk_[m];
			Ti[m] = ik_[m];
		}//H_[m-1] < 1.0e-7
	}//m
	
	if (Nsl == 1)
	{
		if (TV == 5)
		{
			mult(rk_[1],ik_[1], ch[1],ch[2], &u,&v);
			mult(rk_[2],ik_[2], sh[1],sh[2], &u1,&v1);
			u = rk_[1] + rk_[2];
			v = ik_[1] + ik_[2];
		}
		else
			if (TV == 1)
			{
				mult(rk_[1],ik_[1], ch[1],ch[2], &u,&v);
				mult(rk_[2],ik_[2], sh[1],sh[2], &u1,&v1);
				u += u1;
				v += v1;
			}
			else
			{
				mult(rk_[1],ik_[1], sh[1],sh[2], &u,&v);
				mult(rk_[2],ik_[2], ch[1],ch[2], &u1,&v1);
				u += u1;
				v += v1;
			}
		}//Nsl == 1
		else
			if (Nsl == 2)
			{
				u = rk_[1] + rk_[3];
				v = ik_[1] + ik_[3];
				mult(u,v, rk_[2],ik_[2], &u,&v);
				mult(u,v, ch[1],ch[2], &u,&v);

				mult(rk_[1],ik_[1], rk_[3],ik_[3], &u1,&v1);
				mult(rk_[2],ik_[2], rk_[2],ik_[2], &u2,&v2);
				mult(u1 + u2,v1 + v2, sh[1],sh[2], &u1,&v1);
				u += u1;
				v += v1;
		
			}//Nsl == 2
			else
			{
				D[0] = Di[0] = Di[1] =0.0;
				D[1] = 1.0;
				for ( m=1; m <= Nsl-1; m++ )
				{
					mult(D[m-1],Di[m-1], Q[m],Qi[m], &u,&v);
					mult(D[m],Di[m], T[m]+T[m+1],Ti[m]+Ti[m+1], &u1,&v1);
					divide(u1 - u,v1 - v, Q[m+1],Qi[m+1],&u,&v);

					D[m+1] = u;
					Di[m+1] = v;
				}//m
		
				m = Nsl;
				mult(D[m-1],Di[m-1], Q[m],Qi[m], &u,&v);
				mult(D[m],Di[m], T[m]+T[m+1],Ti[m]+Ti[m+1], &u1,&v1);
				u -= u1;
				v -= v1;
			}///Nsl > 2
	
	F_C[1] = u;
	F_C[2] = v;
	
	//printf (" \n	KZ[1]=%5.3e			F_C1=%5.3e",KZ[1],F_C[1] );
	//fprintf( out,"\nKZ[1]=%5.3e			F_C1=%5.3e",KZ[1],F_C[1] );

	//printf (" \n	KZ[1]=%5.3e	KZ[2]=%5.3e		F_C1=%5.3e	F_C2=%5.3e",KZ[1],KZ[2],F_C[1],F_C[2] );
	//fprintf( out,"\n		KZ[1]=%5.3e	KZ[2]=%5.3e		F_C1=%5.3e	F_C2=%5.3e",KZ[1],KZ[2],F_C[1],F_C[2] );

	Ffu = sqrt(u * u + v * v);

	ik++;
	
}/* F */

