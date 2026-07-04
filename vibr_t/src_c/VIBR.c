
#include "VIBRr.h"


int main ( void )
{
	int var, Nvar,q,n,p,m,M_z1,m_[221],m1,P;
	int MZR,mzr,nl;
	double pp, mm, n_m, u, v, kx_, ky_, kz_, Re_eps, Im_eps;
	//double u,v,u1,v1;

	
	/********************
	// Решетка из эллиптических диэлектрических цилиндров 29.08.10

	*******************/

	printf("\n Enter the quantity of variants - ");
	scanf("%d",&Nvar);

	printf("\n Enter the name for input file - ");

	scanf("%s",inf);


	if ((in = fopen(inf,"r")) == NULL)
	{
		printf("\n Error! Input file was not found! \n");
		exit(1);
	}
	printf("\n Enter the name of the file for results - ");

	scanf("%s",outf);

	if ((out = fopen(outf,"a")) == NULL)
	{
		printf("\n The target file is not opened! \n");
		exit(1);
	}

	/**/
	printf("\n ");

	
	printf("\n Will you build the schedules F(frequency)?(Y/N) - ");
	graf2 = yesno();

	if (graf2)
	{
		printf("\n Enter a name of a file for record of the table for construction of diagrams - ");
		scanf("%s",out2f);
		if ((out2 = fopen(out2f,"a")) == NULL)
		{
			printf("\n  The target file is not opened! \n");
			exit(1);
		}
		
	}
	
	

	Pi = 4.0 * atan(1.0);
	Pi_2 = 2.0 * Pi;
	pi2 = 0.5 * Pi;

	SQ_05 = sqrt(0.5);
	Z0 = 120 * Pi;
	var = 0;
	start = clock();
	
	/*
	for (m = 0; m <= 40; m++)
	{
		double l, eps1, eps2;
		l = 1000.0 + 125 * m;
		VO2_cold(l, &eps1, &eps2);
		libmat(20, l, &Re_eps, &Im_eps);
		//printf("\n  l=%6.3lf   EPS_1=%6.3lf    EPS_2=%6.3lf", l, eps1, eps2);
		printf("\n  l=%6.3lf   EPS_1=%6.3lf    EPS_2=%6.3lf", l, eps1 - Re_eps, eps2 - Im_eps);

		VO2_hot(l, &eps1, &eps2);
		libmat(21, l, &Re_eps, &Im_eps);
		//printf("\n  l=%6.3lf   EPS_1=%6.3lf    EPS_2=%6.3lf", l, eps1, eps2);
		printf("\n  l=%6.3lf   EPS_1=%6.3lf    EPS_2=%6.3lf", l, eps1 - Re_eps, eps2 - Im_eps);

		
	}
	fflush(out);
	exit(1);
	*/
/*--------------------------------------------------------------------------*/
	do
	{	
		str = clock();
		var++;
		VVOD(); // input 
		
		AM0( );//new захват памяти
		 		
		E_extern (&kx_,&ky_,&kz_);
		// calc. external downward ElMagn-field
		// A.L. should extend for arbitrary illumination
		
		if (SIM_r > 0)//
			AM_r[1] = AM_r[3] = AM_r[2] = 1;//!!!!!!!!!!!!!!!!!

		for ( N_g = AN_g[1]; N_g <= AN_g[3]; N_g += AN_g[2] ) // LOOP for M_g optimization
		for ( M_z1 = AM_z[1]; M_z1 <= AM_z[3]; M_z1 += AM_z[2] ) // LOOP for M_z1 optimization
		for ( M_fi = AM_fi[1]; M_fi <= AM_fi[3]; M_fi += AM_fi[2] ) // LOOP for M_fi optimization
		for ( M_r_ = AM_r[1]; M_r_ <= AM_r[3]; M_r_ += AM_r[2] )
		{ // LOOP for M_r optimization
			 if (SIM_r >= 0)//
				 M_r_ = M_fi;//
		
		M_r = (N_La == 1) ? M_r_ : M_r_ * N_La + 1;
		M_r_max = M_r;
		M_r_1 = (N_La == 1) ? M_r_ : M_r_ + 1;
		 
		AM4();//new
		
		//ZERO (j,J); // tabulates Zeros of Bessel derivative for 0-8 Bessel  orders for 1-10 first zeroes
		
		//end new 
			
			printf( "\n\n\n\n  N_g=%2d  M_z=%2d      M_fi=%2d  M_r=%2d",N_g,M_z1,M_fi,M_r);
			fprintf( out,"\n\n\n\n  N_g=%2d  M_z=%2d      M_fi=%2d  M_r=%2d",N_g,M_z1,M_fi,M_r);
							
				
			for ( q=0, f = AF[1]; f <= AF[3]; f += AF[2] )
			{ // LOOP over wavelengths f
				q++;
				lam = f;

				if (bum)
				{
					for (n = 1; n <= N_hol; n++)//Ввод диэл. прог. металлических лайнеров
					{
						for (m = 1; m <= N_l[n]; m++)
						{
							if (is_l[n][m])
							{
								libmat(is_l[n][m], f, &Re_eps, &Im_eps);
								EPS_l_r_[n][m] = Re_eps;
								EPS_l_i_[n][m] = Im_eps;

								//printf("\n  EPS_r[m]=%6.3lf   EPS_i[m]=%6.3lf", Re_eps, Im_eps);//
							
							}//is_l[n][m]


						}//m
					}//n
				}//bum

				for (n = 1; n <= N_hol; n++)//Ввод диэл. композитов лайнеров
				{
					for (m = 1; m <= N_l[n]; m++)
					{
						if (comp_l[n][m] > 0.0001)
						{
							COMPOSITE (EPS_l_r_[n][m], EPS_l_i_[n][m], CONCEN_l[n][m], comp_l[n][m], &Re_eps, &Im_eps);

							EPS_l_r[n][m] = Re_eps;
							EPS_l_i[n][m] = Im_eps;
							printf("\n\n\n  Lainer EPS_l_r[n][m]=%6.3lf   EPS_i[m]=%6.3lf\n ",
								EPS_l_r[n][m], EPS_l_i[n][m]);

							COMPOSITE_3(EPS_l_r_[n][m], EPS_l_r_[n][m],
								EPS_l_i_[n][m], EPS_l_i_[n][m],
								CONCEN_l[n][m], 0.0,
								comp_l[n][m], comp_l[n][m],
								&Re_eps, &Im_eps);

							EPS_l_r[n][m] = Re_eps;
							EPS_l_i[n][m] = Im_eps;
							printf("\n  3 Lainer EPS_l_r[n][m]=%6.3lf   EPS_i[m]=%6.3lf\n ",
								EPS_l_r[n][m], EPS_l_i[n][m]);

						}//comp_l[n][m]
						else
						{
							EPS_l_r[n][m] = EPS_l_r_[n][m];
							EPS_l_i[n][m] = EPS_l_i_[n][m];
						}

						//printf("\n  Lainer EPS_l_r[n][m]=%6.3lf   EPS_i[m]=%6.3lf\n ", EPS_l_r[n][m], EPS_l_i[n][m]);//
					}//m
				}//n
				

				if (bum1)
				{
					for (m = 2; m <= Nsl + 1; m++)//Ввод диэл. прог. металлических слоев
					{
						if (is[m])
						{
							libmat(is[m], f, &Re_eps, &Im_eps);

							EPS_r_[m] = Re_eps;
							EPS_i_[m] = Im_eps;
						}//is[n]

						//	printf ( "\n  EPS_r[m]=%6.3lf   EPS_i[m]=%6.3lf",EPS_r[m],EPS_i[m]);//
						//	fprintf ( out,"\n  EPS_r[m]=%6.3lf   EPS_i[m]=%6.3lf",EPS_r[m],EPS_i[m]);
					}//m
				}//bum1


				for (m = 1; m <= Nsl + 1; m++)//Ввод диэл. композитов слоев
				{
					if (comp_[m] > 0.00001)
					{
						COMPOSITE (EPS_r_[m], EPS_i_[m], CONCEN[m], comp_[m], &Re_eps, &Im_eps);

						EPS_r[m] = Re_eps;
						EPS_i[m] = Im_eps;
					}//comp_[n]
					else
					{
						EPS_r[m] = EPS_r_[m];
						EPS_i[m] = EPS_i_[m];
					}

					//	printf ( "\n  EPS_r[m]=%6.3lf   EPS_i[m]=%6.3lf",EPS_r[m],EPS_i[m]);//
					//	fprintf ( out,"\n  EPS_r[m]=%6.3lf   EPS_i[m]=%6.3lf",EPS_r[m],EPS_i[m]);
				}//m


				u = sqrt(fabs(EPS_r[Nsl+1])) / lam;
				Ng_r = (int)(dx * u);
				n =  (int)(dx * u);
				if (Ng_r <  n)
					Ng_r = n;

				for ( M_r = 0, n=1; n <= N_hol; n++ )//new
				{
						MR[n] = (N_l[n] == 1) ? M_r_1 : M_r_ * N_l[n] + 1;//new
						if (M_r < MR[n])//new
							M_r = MR[n];
				}

				// correct number of angular base f-ns for polarization corrected for symmetry caused by illumionation
				for ( P=0; P <= 1; P++ )
				{
					mfi1[P] = M_fi - P;
					mfi2[P] = M_fi - (1 - P);
				
					for ( M_z = MZR = 0, n=1; n <= N_hol; n++ )//new
					{
						m = N_hol - n + 1;
						//n_m = Max(n_[0],n_[m+1]);
						n_m = n_[m+1];

						m_[n] = (M_z1 > 2) ? Max( (int)(H_[m] / f * n_m  * M_z1), 3) : M_z1; // H1/wavelength_corrected*M_z1

						//m_[n] =  2;

						mzr = m_[n] * MR[m];
					
						MZR += mzr;
						M_z += m_[n];

					}//n

				
					m_t[2] = M_z;
					m_b[2] = m_t[2] - m_[N_hol] + 1;
				
					for ( n=3; n <= N_hol+1; n++ )
					{
						m_t[n] = m_b[n-1] - 1;//номер верхней границы слоя
						m_b[n] = m_t[n] - m_[N_hol + 2 - n] + 1;//номер нижней границы слоя
					}
				
					/*				
					ms =  M_z * M_r;
					M_[1] = mfi1[P] * ms; // total number of base f-ns for Ex
					M_[2] = M_[3] = mfi2[P] * ms; // total number of base f-ns for Ey and Ez.
					*/
					M_[1] = mfi1[P] * MZR; // total number of base f-ns for Ex
					M_[2] = M_[3] = mfi2[P] * MZR; // total number of base f-ns for Ey and Ez.

					MS[P] = M_[1] + M_[2] + M_[3];	// total number of base f-ns for all of E
					ms = M_z;
					
				}//P
	
				MS_ = MS[0] + MS[1];

				MS_N = MS_ * N_ell;

				printf( "   M_z_calc=%2d  M_s=%2d",M_z,MS_N);
				fprintf( out,"   M_z_calc=%2d  M_s=%2d",M_z,MS_N);
			
			
			
				printf ("\n\n  lamda =  %5.3lf \n",f);
				fprintf( out,"\n\n  lamda =  %5.3lf \n",f);

				AM5();//new
											
				N_r[0] = 2;
				N_r[ms+1] = N_hol + 1;
				v = 0.0;
				for ( m1 = 0, n=1; n <= N_hol; n++ )
				{
					m = N_hol - n + 1; //N sloya

					// new относительные внещние размеры лайнеров 
					for ( nl=1; nl <= N_l[n]; nl++ )
						r_l[n][nl] = A_t_l[n][nl] / A_t[n];// относительные внещние размеры лайнеров
					
					r_l[n][0] = 0.0;
					
					// new 
				
					h_b[m+1] = v;

					mm = (double)m_[n];
					u = H_[m] / (mm - 1.0);
										
					for ( p=1; p <= m_[n]; p++ )
					{
						m1++;
						N_r[m1] = m + 1;//Номер слоя, в котором расположена точка Zm[m1]
						pp = (double)p;
						Zm[m1] = v + u * (pp - 1.0);
						a_[m1] = 0.5 * (A_t[m] * (pp - 1.0) + A_b[m] * (mm - pp)) / (mm - 1.0);//внещние размеры неоднородности
						b_[m1] = 0.5 * (B_t[m] * (pp - 1.0) + B_b[m] * (mm - pp)) / (mm - 1.0);
						
						for ( nl=1;nl <= N_l[m]; nl++ )//внещние размеры лайнеров
						{
							al_[m1][nl] = 0.5 * (A_t_l[m][nl] * (pp - 1.0) + A_b_l[m][nl] * (mm - pp)) / (mm - 1.0);
							bl_[m1][nl] = 0.5 * (B_t_l[m][nl] * (pp - 1.0) + B_b_l[m][nl] * (mm - pp)) / (mm - 1.0);
						
						}
							
					}//p
					v = Zm[m1];
				}//n

				hm1[1] = hm1[ms+1] = 0.0;
				for ( m=2; m <= ms; m++ )
				{
					hm[m] = Zm[m] - Zm[m-1];
					hm1[m] = (hm[m] > 1.0e-7) ? 1.0 / hm[m] : 0.0;
					
				}
				
				ZERO_B();//new
				
				N_corn = 2;
				Alap = AllocC(MS_N,MS_N); // alloc memory for matrix inversion using MKL
				Blap = AllocC(N_corn,MS_N); // alloc memory for right side for matrix inversion using MKL

				k = k0 = 2.0 * Pi / f;
				k2 = k * k;
								
				for ( m=1; m <= Nsl + 1; m++ )
				{
					k2_r[m] = k2 * EPS_r[m];
					k2_i[m] = k2 * EPS_i[m];
		
					k_r[m] = k * n_[m];
					k_i[m] = -k * k_[m]; // "-" minus sign to aliign with exp(iwt) 
				}

				kz = k * kz_;
				ky = k * ky_;
				kx = k * kx_;	// due to illum.azimuth along Y-axis		
				
				EXTERN(E,H);  // calc. El.Magn field accounting for refl.&propagation through 5 solid layers

				//if (H_[0] < 0.5)
				{
					printf ("\n  p  q    Rpq    Tpq         delta P");                                                      
					fprintf( out,"\n  p  q     Rpq         Tpq            delta P");
			
				}
				/*
				else
				{
					printf ("\np  q     Rpq         Tpq             Psum        Psum * Tp       delta P");                                                      
					fprintf( out,"\n  q     Rpq         Tpq             Psum        Psum * Tp       delta P");                                                      
				}
				*/
		

				F_df(q); //new освобождение памяти

				FM5( );
						
			}/* f */ // LOOP over wavelengths f
			
			FM4( );//new	
					
			if (graf2)
				fprintf( out2,"\n     fffffffffff");


			if (graf1)
				fprintf( out1,"\n     fffffffffff");
			
	 // LOOP for M_r optimization	
	}/* N_g */ 
		
	FM1( ); 
	FM2( );
	FM3( );
	FM0( );//new освобождение памяти	
	fflush( out );
	
	/*tim ( str);*/
	tim ( str);
	printf ( "\n\n\n__________________________________________________________________");
	fprintf ( out,"\n\n\n__________________________________________________________________");

	} while( var <  Nvar );

		
	if ( Nvar > 1 )
		tim ( start);
	return ( 0 );
}   /*  end main  */


double jn_m(int n,double x)
{
	double dn;
	
	dn = (double)n;
	return ( dn * log(0.5 * x / dn) < -68.0 ? 0.0 : jn(n,x) );
}//jn_m

void COMPOSITE(double R_eps_m, double I_eps_m, double C, int n_mat, double *R_eps_c, double *I_eps_c)
{
	double Re_eps, Im_eps, rh,ih,rp, ip, u, v;

	libmat(n_mat, f, &Re_eps, &Im_eps);

	u = Re_eps - R_eps_m;
	v = Im_eps - I_eps_m;
	rh = Re_eps + 2.0 * R_eps_m;
	ih = Im_eps + 2.0 * I_eps_m;
	divide_(u, v, rh, ih, &rp, &ip);

	u = 1.0 + 2.0 * C * rp;
	v = 2.0 * C * ip;
	rh = 1.0 - C * rp;
	ih = - C * ip;
	divide_(u, v, rh, ih, &rp, &ip);

	mult_(R_eps_m, I_eps_m, rp,ip, &u, &v);
	*R_eps_c = u;
	*I_eps_c = v;

}//COMPOSITE



void COMPOSITE_3(double R_eps_m_1, double R_eps_m_2, 
	double I_eps_m_1, double I_eps_m_2, 
	double C1, double C2,
	int n_mat_1, int n_mat_2, 
	double* R_eps_c, double *I_eps_c)
{
	double Re_eps_1, Im_eps_1, rh1, ih1, rp1, ip1, u1, v1;
	double Re_eps_2, Im_eps_2, rh2, ih2, rp2, ip2, u2, v2;
	double Re_sum, Im_sum;

	libmat(n_mat_1, f, &Re_eps_1, &Im_eps_1);
	libmat(n_mat_2, f, &Re_eps_2, &Im_eps_2);

	u1 = Re_eps_1 - R_eps_m_1;	//eps_p_1-eps_m Re
	v1 = Im_eps_1 - I_eps_m_1;	//eps_p_1-eps_m Im
	rh1 = Re_eps_1 + 2.0 * R_eps_m_1;	//eps_p_1+2*eps_m Re
	ih1 = Im_eps_1 + 2.0 * I_eps_m_1;	//eps_p_1+2*eps_m Im
	divide_(u1, v1, rh1, ih1, &rp1, &ip1);

	u2 = Re_eps_2 - R_eps_m_2;	//eps_p_1-eps_m Re
	v2 = Im_eps_2 - I_eps_m_2;	//eps_p_1-eps_m Im
	rh2 = Re_eps_2 + 2.0 * R_eps_m_2;	//eps_p_1+2*eps_m Re
	ih2 = Im_eps_2 + 2.0 * I_eps_m_2;	//eps_p_1+2*eps_m Im
	divide_(u2, v2, rh2, ih2, &rp2, &ip2);

	Re_sum =C1 * rp1 + C2 * rp2;
	Im_sum = C1 * ip1 + C2 * ip2;

	u1 = 1 + 2 * Re_sum;
	v1 = 2 * Re_sum;
	rh1 = 1.0 - Re_sum;
	ih1 = -1.0 * Im_sum;
	divide_(u1, v1, rh1, ih1, &rp1, &ip1);
	mult_(R_eps_m_1, I_eps_m_1, rp1, ip1, &u1, &v1);

	*R_eps_c = u1;
	*I_eps_c = v1;

	/*
	u1 = 1.0 + 2.0 * C1 * rp1;
	v1 = 2.0 * C1 * ip1;
	rh1 = 1.0 - C1 * rp1;
	ih1 = -C1 * ip1;
	divide_(u1, v1, rh1, ih1, &rp1, &ip1);

	u2 = 1.0 + 2.0 * C2 * rp2;
	v2 = 2.0 * C2 * ip2;
	rh2 = 1.0 - C2 * rp2;
	ih2 = -C2 * ip2;
	divide_(u2, v2, rh2, ih2, &rp2, &ip2);

	mult_(R_eps_m_1, I_eps_m_1, rp1, ip1, &u1, &v1);
	mult_(R_eps_m_2, I_eps_m_2, rp2, ip2, &u2, &v2);
	*R_eps_c = u1+u2;
	*I_eps_c = v1+u2;
	*/
}//COMPOSITE


void COMPOSITE_me(double R_eps, double I_eps, //эпс матрицы
	double C,// концентрация наночастиц
	int n_mat, // номер диэлектрика наночастиц
	double *R_eps_c, double *I_eps_c// эффективная эпс
	)
{
	double Re_eps, Im_eps, rh, ih, u, v;

	libmat(n_mat, f, &Re_eps, &Im_eps);

	u = Re_eps - R_eps;
	v = Im_eps - I_eps;

	mult_(3.0 * R_eps, 3.0 * I_eps, u, v, &u, &v);

	rh = Re_eps + 2.0 * R_eps;
	ih = Im_eps + 2.0 * I_eps;

	divide_(u, v, rh, ih, &rh, &ih);

	*R_eps_c = 1.0 + (1.0 - C) * (R_eps - 1.0) + C * rh;
	*I_eps_c = (1.0 - C) * I_eps + C * ih;

}//COMPOSITE