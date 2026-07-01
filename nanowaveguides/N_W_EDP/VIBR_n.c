
#include "VIBRr.h"

void EPS_MET(int is, double *EPS_r, double *EPS_i);

int main ( void )
{
	int var, Nvar,q,n,p,m,M_z1,M_r1,m_z[21],m1,mu,nu,N_g1,mh,nn,mm;
	double pp,m_,n_m,u,v,kx_,ky_,kz_,y1,y2;
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

		if (Pol)
		{
			i_c = 2;
			i_f = 3;
		}
		else
			i_c = i_f = 1;
			
		E = ( double****)amat_4d( 2,N_hol_m+1,1,2,1,3, 1,4, sizeof( double) );
		H = ( double****)amat_4d( 2,N_hol_m+1,1,2,1,3, 1,4, sizeof( double) );
		Ep = ( double****)amat_4d( 2,N_hol_m+1,1,2,1,3, 1,4, sizeof( double) );
		Hp = ( double****)amat_4d( 2,N_hol_m+1,1,2,1,3, 1,4, sizeof( double) );
		

		j = ( double**)amat_2d( 0,8, 1,10, sizeof( double) );
		J = ( double**)amat_2d( 0,8, 1,10, sizeof( double) );
						
		E_extern (&kx_,&ky_,&kz_);
		// calc. external downward ElMagn-field
		// A.L. should extend for arbitrary illumination
		
	
		for ( N_g1 = AN_g[1]; N_g1 <= AN_g[3]; N_g1 += AN_g[2] ) // LOOP for M_g optimization
		for ( M_z1 = AM_z[1]; M_z1 <= AM_z[3]; M_z1 += AM_z[2] ) // LOOP for M_z1 optimization
		for ( M_r1 = AM_r[1]; M_r1 <= AM_r[3]; M_r1 += AM_r[2] )
		{ // LOOP for M_r optimization
			M_fi = 1;
			

			printf( "\n\n  N_g=%2d  M_z=%2d    M_r=%2d",N_g1,M_z1,M_r1);
			fprintf( out,"\n\n  N_g=%2d  M_z=%2d    M_r=%2d",N_g1,M_z1,M_r1);
							
				
			for ( q=0, f = AF[1]; f <= AF[3]; f += AF[2] )
			{ // LOOP over wavelengths f
				q++;
				lam = f;
				k0 = 2.0 * Pi / f;

			
				N_g = (int)(.2 * sqrt(fabs(EPS_r[0])) * k0 * dy) + N_g1;
			
				for ( m=0; m <= 5; m++ )//????????????????????????????????????
				{
					if (!is[m])
						continue;

					// for metals START
					EPS_MET(is[m], &u,&v);
					EPS_r[m] = u;
					EPS_i[m] = v;
					
					/*
					if ( is[m] == 2 )
					{
						EPS_r[m] = Ag(f,0);
						EPS_i[m] = Ag(f,1);
					}
					else 
						if ( is[m] == 3 )
						{
							EPS_r[m] = Au(f,0);
							EPS_i[m] = Au(f,1);
						}
						else
							if ( is[m] == 1 )
							{
								EPS_r[m] = Cu(f,0);
								EPS_i[m] = Cu(f,1);
							}
					*/

				//	printf ( "\n  EPS_r[m]=%6.3lf   EPS_i[m]=%6.3lf",EPS_r[m],EPS_i[m]);//
				//	fprintf ( out,"\n  EPS_r[m]=%6.3lf   EPS_i[m]=%6.3lf",EPS_r[m],EPS_i[m]);//
					// for metals END
				}

				if (is_in)
				{ // for metals 
					EPS_MET(is_in, &u,&v);
					EPS_r_in = u;
					EPS_i_in = v;
				}
					
				// correct number of angular base f-ns for polarization corrected for symmetry caused by illumionation
				for ( ms = M_z = M_r = 0, n=1; n <= N_hol; n++ )
				{
					mm = N_hol - n + 1;

					nn = N_discon[n]; // номер слоя с n -ой неоднородностью
					
					//n - нумерация неоднородностей идет сверху вниз, например, 2, 4, 6.
					//нумерация сплайнов по z идет снизу вверх.
					
					n_m = Max(n_[0],n_[nn]);
					m_z[mm] =  Max( (int)(H_[nn-1] / f * n_m  * M_z1), 3); // число узлов по z в слое с n-ой неоднородностью
					M_z += m_z[mm];//  общее число узлов по z
					
				//	printf (" \n n= %2d  mz= %2d  ",mm,m_z[mm]);
				//	fprintf( out,"\n n= %2d  mz= %2d  ",mm,m_z[mm]);
				
					N_di[nn] = n;
					m_r[mm] = 0;// общее число узлов по r в слое с n-ой неоднородностью
					for ( M_r_max = nu = 0,m=1;m <= N_l[mm]; m++ )
					{
						Dmax = 0.5 * (B_t[mm][m] + B_b[mm][m]);
						m_r_[mm][m] = (int)(Dmax * k0 * .25 * sqrt(fabs( EPS_l_r[mm][m]))) + M_r1;// число узлов по r в m-о зубе,
																														  //лежащем в слое с mm-ой неоднородностью
						m_r[mm] += m_r_[mm][m];
					
						if (m_r_[mm][m] > M_r_max)
							M_r_max = m_r_[mm][m]; 

						/*
						for ( mu=1;mu <= m_r[mm]; mu++ )//UBRAT
						{
							nu++;
							NRl[nn][nu] = m;//по номеру БФ по r (nu) и номеру слоя с неоднородностью (nn) восстанавливает номер зуба
						}
						*/

					}//m
					
					//printf( "   n=%2d  mm=%2d =%2d    N_l[mm]=%2d m_z[mm]=%2d  ",N_l[mm],m_z[mm]);
					//fprintf( out,"   M_z_calc=%2d \n ",M_z);
					M_r += m_r[mm];
					
					m_s[mm] = m_r[mm] * m_z[mm];//число БФ в слое с n-ой неоднородностью

					ms += m_s[mm];//общее число БФ
				}//n

				
				MS_ = (Pol) ? 2 * ms : ms;
					
				MS__ = ms;
				ms = M_z;//????????????????????
									
				
				printf ("\n\n  lamda =  %5.3lf ",f);
				fprintf( out,"\n\n  lamda =  %5.3lf ",f);

				printf( "   M_z_calc=%2d \n ",M_z);
				fprintf( out,"   M_z_calc=%2d \n ",M_z);
							
				// memory allocation for 1D-array
				ASR_ = ( double**)amat_2d(0,1,1,MS_, sizeof( double) );// third party code (source belongs to A.L. colleague)
				ASI_ = ( double**)amat_2d(0,1,1,MS_, sizeof( double) );
				Zm = ( double*)amat_1d( 0,ms+1, sizeof( double) );
				hm = ( double*)amat_1d( 1,ms+1, sizeof( double) );
				hm1 = ( double*)amat_1d( 1,ms+1, sizeof( double) );
				dm = ( double*)amat_1d( 1,ms, sizeof( double) );
				
				b_ = ( double**)amat_2d( 1,ms, 1,N_max,sizeof( double) );
				Y_centr = ( double**)amat_2d( 1,ms, 1,N_max,sizeof( double) );
				MZ_C = ( int*)amat_1d( 1,ms, sizeof( int ) );
				MZ_F = ( int*)amat_1d( 1,ms, sizeof( int ) );
				
						
				N_r = ( int*)amat_1d( 0,ms+1, sizeof( int) );
				// memory allocation for 4D-array
				No = ( int***)amat_3d( 1,3, 1,M_r, 1,ms, sizeof( int) );//??????????????????????
				No2 = ( int***)amat_3d( 1,3, 1,M_r, 1,ms, sizeof( int) );
				
            	BAS = ( double***)amat_3d( 0,N_g, 1,M_z, 1,M_r, sizeof( double) );
				
				if (P_S)
				{
					BAS_z = ( double***)amat_3d( 0,N_g, 1,M_z, 1,M_r, sizeof( double) );
					EVEN_B = ( int**)amat_2d(  1,M_z, 1,M_r, sizeof( int) );
				}


				{
					//Блок вычисления узлов сплайнов по координате z
					//Начало координат z=0 - нижняя граница слоя с неоднородностью с № N_hol
					double B[21],del;
					//B[n] - верхняя граница -го слоя. Нумерация снизу вверх. Начало нумерации от z=0.
					
					int Ncom;//

					Ncom = N_discon[N_hol] - 1;

					B[0] = 0;
					for ( p=0; p < Ncom; p++ )
						B[p+1] = B[p] + H_[Ncom - p];//proveril

					Ncom++;

					N_r[0] = 2;
					N_r[ms+1] = N_hol + 1;
					v = 0.0;
					for ( m1 = 0, n=1; n <= N_hol; n++ )
					{
						mh = N_hol - n + 1; // нумерация слоев от плоскости z=0 - нижней грани нижнего слоя нижней неоднородности. Её номер N_hol
						m = N_discon[mh]; //N sloya  с неоднородностью
				
						//h_b[m+1] = v;
						
						//nu =  N_discon[mh] - 1;
						nu = Ncom - m;
						h_b[m] =  B[nu];

						m_ = (double)m_z[n];
					//	u = H_[nu] / (mm - 1.0);
						del = (B[nu+1] - B[nu]) / (m_ - 1.0);;

						for ( p=1; p <= m_z[n]; p++ )
						{
							m1++;
							N_r[m1] = m;//Номер слоя, в котором расположена точка Zm[m1]
							pp = (double)p;
						
							Zm[m1] = B[nu] + del * (pp - 1.0);
							
							for ( mu=1;mu <= N_l[mh]; mu++ )
							{
								y1 = (Ytc[mh][mu] * (pp - 1.0) + Ybc[mh][mu] * (m_ - pp)) / (m_ - 1.0);
								y2 = (Ytf[mh][mu] * (pp - 1.0) + Ybf[mh][mu] * (m_ - pp)) / (m_ - 1.0);
								b_[m1][mu] = 0.5 * (y2 - y1);
								Y_centr[m1][mu] = 0.5 * (y2 + y1);
	
							}//m
						}//p
						v = Zm[m1];
					}//n
					

					hm1[1] = hm1[ms+1] = 0.0;
					for ( m=2; m <= ms; m++ )
					{
						hm[m] = Zm[m] - Zm[m-1];
						hm1[m] = (hm[m] > 1.0e-7) ? 1.0 / hm[m] : 0.0;
					
					}

					//m = 2;
					m = N_discon[1];
					m_t[m] = M_z;
					m_b[m] = m_t[m] - m_z[N_hol] + 1;
				
					mu = N_hol;
				
					for ( n=3; n <= N_hol+1; n++ )
					{
						//m = n;
						m = N_discon[n-1];
						m1 = N_discon[n-2];
						mu = N_hol + 2 - n;
						m_t[m] = m_b[m1] - 1;//номер верхней границы слоя
						m_b[m] = m_t[m] - m_z[mu] + 1;//номер нижней границы слоя
					}
				
			
				}//Блок вычисления узлов сплайнов по координате z. END
				
				
				N_corn = 1;
				Alap = AllocC(MS_,MS_); // alloc memory for matrix inversion using MKL
				Blap = AllocC(N_corn,MS_); // alloc memory for right side for matrix inversion using MKL

				k = 2.0 * Pi / f;
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

				printf ("\n  p  q    Rpq     Tpq       delta P");                                                      
				fprintf( out,"\n  p  q    Rpq     Tpq       delta P");
				
                BASIS_L();

				F_df(q); // !!! core f-n

				DeAllocC(Alap);
				DeAllocC(Blap);
			
				fmat_1d( (void*)MZ_F,  1, sizeof( int) );
				fmat_1d( (void*)MZ_C,  1, sizeof( int) );
				fmat_1d( (void*)N_r,  0, sizeof( int) );
				fmat_2d( (void**)b_, 1,ms, 1, sizeof( double) );
				fmat_2d( (void**)Y_centr, 1,ms, 1, sizeof( double) );
				fmat_2d( (void*)ASI_, 0,1, 1, sizeof( double) );//
				fmat_2d( (void*)ASR_, 0,1, 1, sizeof( double) );//
				fmat_1d( (void*)Zm,  0, sizeof( double) );
				fmat_1d( (void*)hm,  1, sizeof( double) );
				fmat_1d( (void*)hm1,  1, sizeof( double) );
				fmat_1d( (void*)dm,  1, sizeof( double) );
				fmat_3d( (void***)No,1,3, 1,M_r,   1, sizeof( int) );
				fmat_3d( (void***)No2,1,3, 1,M_r,   1, sizeof( int) );
	
           		fmat_3d( (void***)BAS,0,N_g, 1,M_z, 1, sizeof( double) );
				
				if (P_S)
				{
					fmat_3d( (void***)BAS_z,0,N_g, 1,M_z, 1, sizeof( double) );
					fmat_2d( (void**)EVEN_B, 1,M_z, 1, sizeof( int) );
				}
 
								
			}/* f */ // LOOP over wavelengths f
			
			
			if (graf2)
				fprintf( out2,"\n     fffffffffff");


			if (graf1)
				fprintf( out1,"\n     fffffffffff");
			
	 // LOOP for M_r optimization	
	}/* N_g */ 
		
		
		
		fflush( out );

		/*tim ( str);*/

		tim ( str);
		printf ( "\n\n\n__________________________________________________________________");
		fprintf ( out,"\n\n\n__________________________________________________________________");

	} while( var <  Nvar );

	fmat_4d( (void****)E, 2,N_hol_m+1,1,2,1,3, 1, sizeof( double) );
	fmat_4d( (void****)H, 2,N_hol_m+1,1,2,1,3, 1, sizeof( double) );
	fmat_4d( (void****)Ep, 2,N_hol_m+1,1,2,1,3, 1, sizeof( double) );
	fmat_4d( (void****)Hp, 2,N_hol_m+1,1,2,1,3, 1, sizeof( double) );
	
			
	fmat_2d( (void**)j, 0,8, 1, sizeof( double) );
	fmat_2d( (void**)J, 0,8, 1, sizeof( double) );
	
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

void EPS_MET(int is, double *EPS_r, double *EPS_i)
{
	if ( is == 2 )
	{
		*EPS_r = Ag(f,0);
		*EPS_i = Ag(f,1);
	}
	else 
		if ( is == 3 )
		{
			*EPS_r = Au(f,0);
			*EPS_i = Au(f,1);
		}
		else
			if ( is == 1 )
			{
				*EPS_r = Cu(f,0);
				*EPS_i = Cu(f,1);
			}
}