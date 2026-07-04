
#include "VIBRr.h"

void EXTERN(double ****E, double ****H)//внешнее поле = падающее поле + отраженное от многослойного диэлектрика
{ // calc. El.Magn field accounting for refl.&propagation through 5 solid layers
	double u,v,u1,v1,u2,v2,x,y,t[111],it[111],fi[111],ifi[111],C[111],iC[111],u_,v_,
		kap_i[111],kap_r[111],kap1_i[111],kap1_r[111],q[111],iq[111],z,T,iT,
		rE2,iE2,rE3,iE3,Z1,s[111],si[111],EPS_ii[111],y_[111],pci[111],ex;
	int N_d,n,i,j;

	//N_d - nomer sloya, v kotorom neodnorodnost.


	/* RESULT
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


	rX[1],iX[1] - R complex reflection coeff in ambient 
	rX[2],iX[2] - T complex transmission coeff. in substrate
	rX[3],iX[3] - A complex.coeff of field within layer #2 without insertion
	rX[4],iX[4] - B complex.coeff of field within layer #2 without insertion
	*/
	for ( N_d = 2; N_d <= N_hol + 1; N_d++ )
		for ( n=1; n <= 2; n++ )
			for ( i=1; i <= 3; i++ )
				for ( j=1; j <= 4; j++ )
					E[N_d][n][i][j] = H[N_d][n][i][j] = Ep[N_d][n][i][j] = Hp[N_d][n][i][j] = 0.0;
		


	
	for ( n=1; n <= Nsl+1; n++ )
		EPS_ii[n] = EPS_i[n] - 1.0e-19;
	
	for ( Pol=0; Pol <= 1; Pol++ )
	{
	for ( n=1; n <= Nsl+1; n++ )// вычисления капа calculating the cap
	{
		x = k2_r[n] - ky * ky - kx * kx;
		y = k2_i[n];
		
		C_sqrt(x,y,&u,&v);
		kap_r[n] = u;
		kap_i[n] = v;
		
		if (Pol)
		{
			divide_(kap_r[n],kap_i[n], EPS_r[n],EPS_ii[n], &u,&v); // divide 2 complex numbers
			kap1_r[n] = u;
			kap1_i[n] = v;
		}
		else
		{
			kap1_i[n] = kap_i[n];
			kap1_r[n] = kap_r[n];
		}
	}//n

	// вачисления коэффициентов A, B, R, T calculate the coefficients
// _______________________________________________________
	for ( n=2; n <= Nsl; n++ )
	{
		x = H_[n-1] * kap_r[n];
		y = y_[n] = H_[n-1] * kap_i[n];
		
		//C_exp(y,-x,&u1,&v1);
		u1 = cos(x);
		v1 = -sin(x);
		
		x *= 2.0;
		y *= 2.0;
		C_exp(y,-x,&u2,&v2);

		divide_(1.0 + u2,v2, 1.0 - u2,-v2, &u,&v);
		mult_(-kap1_i[n],kap1_r[n], u,v, &t[n],&it[n]);//t

		divide_(-2.0 * v1,2.0 * u1, 1.0 - u2,-v2, &u,&v);
		ex = exp(y);
		
		s[n] = u * ex;
		si[n] = v * ex;

		mult_(kap1_r[n],kap1_i[n], u,v, &q[n],&iq[n]);//q
	}

	t[1] = -kap1_i[1];
	it[1] = kap1_r[1];
	t[Nsl+1] = -kap1_i[Nsl+1];
	it[Nsl+1] = kap1_r[Nsl+1];

	q[1] = iq[1] = q[Nsl+1] = iq[Nsl+1] = 0.0;

	for ( n=1; n <= Nsl; n++ )
	{
		fi[n] = t[n] + t[n+1];
		ifi[n] = it[n] + it[n+1];
	}

	C[Nsl+1] = 1.0;

	y_[Nsl+1] = pci[Nsl+1] = pci[Nsl+2] = iC[Nsl+1] = iC[Nsl+2] = C[Nsl+2] = 0.0;

	for ( n=Nsl; n >= 2; n-- )
		pci[n] = pci[n+1] + y_[n];

	for ( n=Nsl; n >= 2; n-- )
	{
		mult_(C[n+1],iC[n+1], fi[n],ifi[n], &u,&v);
		mult_(C[n+2],iC[n+2], q[n+1],iq[n+1], &u1,&v1);
		ex = exp(2.0 * y_[n+1]);
		divide_(u - u1 * ex,v - v1 * ex, q[n],iq[n], &u1,&v1);
		C[n] = u1;
		iC[n] = v1;
	}//n

	mult_(C[2],iC[2], fi[1],ifi[1], &u,&v);
	mult_(C[3],iC[3], q[2],iq[2], &u1,&v1);
	ex = exp(2.0 * y_[2]);
	divide_(2.0 * t[1],2.0 * it[1],u - u1 * ex,v - v1 * ex,  &T,&iT);

	for ( n=2; n <= Nsl+1; n++ )
	{
		mult_(C[n],iC[n], T,iT, &u,&v);
		ex = exp(pci[2] - pci[n]);
		C[n] = u * ex;//A
		iC[n] = v * ex;
	}
	
	ex = exp(pci[2]);
	T *= ex;
	iT *= ex;
	//Вычисление T END calculate______________________________________________________
	
	rX[2] = T;
	iX[2] = iT;
	
	rX[1] = C[2] - 1.0;// коэффициент R factor
	iX[1] = iC[2];
		
	u_ = kap1_r[1] * 2.0 * t_;
	v_ = sin(u_);
	u_ = cos(u_);
	
	for ( j=1; j <= 2; j++ )
		mult_( rX[j],iX[j], u_,v_, &rX[j],&iX[j]);
		
		
	for ( N_d = 2; N_d <= N_hol + 1; N_d++ )
	{
		rX[3] = C[N_d];	
		iX[3] = iC[N_d];

		rX[4] = C[N_d+1];	
		iX[4] = iC[N_d+1];

				
		for ( j=3; j <= 4; j++ )
		{
			mult_( rX[j],iX[j], u_,v_, &rX[j],&iX[j]);
			mult_( rX[j],iX[j], s[N_d],si[N_d], &rX[j],&iX[j]);
		}
		/*	n=1;
		for ( j=1; j <= 4; j++ )
				{
					printf ("\n	  %1d	   %1d	 %1d	 %1d	 \n  %7.5e       %7.5e",
								Pol,N_d,j,n,rX[j], iX[j] );
							fprintf( out,"\n	  %1d	   %1d	 %1d	 %1d	 \n  %7.5e       %7.5e",
								Pol,N_d,j,n,rX[j], iX[j] );
				}
				*/
		// End вачисления коэффициентов calculate A, B, R, T____________________________________________________
			/* RESULT
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


			rX[1],iX[1] - R complex reflection coeff in ambient 
			rX[2],iX[2] - T complex transmission coeff. in substrate
			rX[3],iX[3] - A complex.coeff of field within layer #2 without insertion
			rX[4],iX[4] - B complex.coeff of field within layer #2 without insertion
			*/	
	
		Z1 = Z0 / sqrt(EPS_r[1]);//???????????????????????????????
		//Z1 = Z0;

		if (Pol)//p-поляризация polarization
		{
			z = - sitet * EPS_r[1];
			for ( j=1; j <= 4; j++ )
			{
				Hp[N_d][1][1][j] = rX[j] * sifi / Z1;//Hx!!
				Hp[N_d][2][1][j] = iX[j] * sifi / Z1;//!!

				Hp[N_d][1][2][j] = -rX[j] * cofi / Z1;//Hx//!!
				Hp[N_d][2][2][j] = -iX[j] * cofi / Z1;//!!

				n = (j == 1) ? 1 : (j == 2) ? Nsl+1 : N_d;	
				divide_(rX[j] * z,iX[j] * z, EPS_r[n],EPS_i[n],  &rE2,&iE2);
			
				Ep[N_d][1][3][j] = rE2;//Ez
				Ep[N_d][2][3][j] = iE2;
				
			}//j

			z =  sqrt(EPS_r[1]) / k0;

			Ep[N_d][1][1][1] = -rX[1] * cotet * cofi;//!!
			Ep[N_d][2][1][1] = -iX[1] * cotet * cofi;//!!
				
			Ep[N_d][1][2][1] = -rX[1] * cotet * sifi;//!!
			Ep[N_d][2][2][1] = -iX[1] * cotet * sifi;//!!


			mult_(rX[2],iX[2], kap1_r[Nsl+1],kap1_i[Nsl+1],  &rE2,&iE2);
			
			Ep[N_d][1][1][2] = rE2 * cofi * z;//Ex prosh??????
			Ep[N_d][2][1][2] = iE2 * cofi * z;
		
			Ep[N_d][1][2][2] = rE2 * sifi * z;//Ey prosh//!!
			Ep[N_d][2][2][2] = iE2 * sifi * z;
		
			for ( j=3; j <= 4; j++ )
			{
				mult_(rX[j],iX[j], kap1_r[N_d],kap1_i[N_d],  &rE2,&iE2);
				rE3 = iE2 * z;
				iE3 = - rE2 * z;
				
				Ep[N_d][1][1][j] = rE3 * cofi;//Ex sloy N_d 
				Ep[N_d][2][1][j] = iE3 * cofi;

				Ep[N_d][1][2][j] = rE3 * sifi;//Ey sloy N_d//!!
				Ep[N_d][2][2][j] = iE3 * sifi;
			}
			
		}//Pol
		else
		//s-поляризация polarization
		{
			//z = sitet * n_[1] / Z0;
			for ( j=1; j <= 4; j++ )
			{
				E[N_d][1][1][j] = rX[j] * sifi;
				E[N_d][2][1][j] = iX[j] * sifi;

				E[N_d][1][2][j] = -rX[j] * cofi;
				E[N_d][2][2][j] = -iX[j] * cofi;

				// Z компонента магнитного поля в слоях не нужна
				//H[N_d][1][3][j] = rX[j] * z;
				//H[N_d][2][3][j] = iX[j] * z;
	
			}//j

			z = -1.0 / (k0 * Z0);

			mult_(rX[1] * z,iX[1] * z, kap1_r[1],kap1_i[1],  &rE2,&iE2);

			H[N_d][1][1][1] = - rE2 * cofi;
			H[N_d][2][1][1] = - iE2 * cofi;
			
			H[N_d][1][2][1] = - rE2 * sifi;
			H[N_d][2][2][1] = - iE2 * sifi;
			
			mult_(rX[2] * z,iX[2] * z, kap1_r[Nsl+1],kap1_i[Nsl+1],  &rE2,&iE2);
			H[N_d][1][1][2] = rE2 * cofi;
			H[N_d][2][1][2] = iE2 * cofi;

			H[N_d][1][2][2] = rE2 * sifi;
			H[N_d][2][2][2] = iE2 * sifi;

			/*	магнитное поле в слоях не нужно
			for ( j=3; j <= 4; j++ )
			{
				mult_(rX[j],iX[j], kap1_r[N_d],kap1_i[N_d],  &rE2,&iE2);
				rE3 = iE2;
				iE3 = - rE2;
				H[N_d][1][2][j] = rE3 * z;
				H[N_d][2][2][j] = iE3 * z;
			}
			*/
				
		}//!P
	}//N_d
	}//Pol

}/* EXTERN */


