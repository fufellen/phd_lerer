
#include "VIBRr.h"

void EPS_MET(int is, double *EPS_r, double *EPS_i);


int main ( void )
{
	int var, Nvar,q,n,p,m,M_z1,m_[221],m1,P;
	int MZR,mzr,nl;
	double pp,mm,n_m,u,v,kx_,ky_,kz_;
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

