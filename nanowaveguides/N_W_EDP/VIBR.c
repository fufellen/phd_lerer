
#include "VIBRr.h"

void EPS_MET(int is, double *EPS_r, double *EPS_i);

int main ( void )
{
	int var, Nvar,m,n,i;
		//,m,M_z1,M_r1,m_z[151],m1,mu,nu,N_g1,mh,nn,mm,N_c;
	double v,u,*bet_r,*bet_i,eps_ef_r[21][21],eps_ef_i[21][21];
	
	
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
		
		bet_r = ( double*)amat_1d( 1,Nkor, sizeof( double) );
		bet_i = ( double*)amat_1d( 1,Nkor, sizeof( double) );
		
		for ( f = AF[1]; f <= AF[3]; f += AF[2] )
		{ 
			k = Pi_2 / f;
			k2 = k * k;

			EPS_min = 100000.0;
			EPS_MAX = 0.0;

			
			for ( n=1; n <= Ny; n++ )//Ввод диэл. прог. металлических слоев
				for ( m=1; m <= Nx_; m++ )
				{ 
					if (bum1)
					{
						if (is[n][m])
							{
								if ( is[n][m] == 2 )
								{
									//eps_r[n][m] = Ag(f,0);
									//eps_i[n][m] = Ag(f,1);
									Ag_(f,&eps_r[n][m],&eps_i[n][m]);
								}
								else 
									if ( is[n][m] == 3 )
									{
										//eps_r[n][m] = Au(f,0);
										//eps_i[n][m] = Au(f,1);
										Au_(f,&eps_r[n][m],&eps_i[n][m]);
										//eps_r[n][m] = -28.28;
										//eps_i[n][m] = - 1.557;
									}
									else
										if ( is[n][m] == 1 )
										{
											eps_r[n][m] = Cu(f,0);
											eps_i[n][m] = Cu(f,1);
										}

								eps_r_[n][m] = eps_r[n][m];
								eps_i_[n][m] = eps_i[n][m];
					
							}//is[n]
					}//bum1

					if (EPS_MAX < fabs(eps_r[n][m]))
						EPS_MAX = fabs(eps_r[n][m]);

					if (EPS_min > fabs(eps_r[n][m]))
						EPS_min = fabs(eps_r[n][m]);
				}//n
			

			
				for ( n=1; n <= nsub; n++ )//Ввод диэл. прог. металлических слоев подложек
				{ 
					if (bum2)
					{
						if (is_sub[n])
							{
								if (is_sub[n] == 2 )
								{
									//eps_sub_r[n] = Ag(f,0);
									//eps_sub_i[n] = Ag(f,1);
									Ag_(f,&eps_sub_r[n],&eps_sub_i[n]);
								}
								else 
									if ( is_sub[n] == 3 )
									{
										//eps_sub_r[n] = Au(f,0);
										//eps_sub_i[n] = Au(f,1);
										Au_(f,&eps_sub_r[n],&eps_sub_i[n]);
									}
									else
										if ( is_sub[n] == 1 )
										{
											eps_sub_r[n] = Cu(f,0);
											eps_sub_i[n] = Cu(f,1);
										}
													
							}//is[n]
					}//bum2
					
					if (EPS_MAX < fabs(eps_sub_r[n]))
							EPS_MAX = fabs(eps_sub_r[n]);

						if (EPS_min > fabs(eps_sub_r[n]))
							EPS_min = fabs(eps_sub_r[n]);

				} //n 
			
			
						
			printf ("\n\n  lamda =  %5.3lf ",f);
			fprintf( out,"\n\n  lamda =  %5.3lf ",f);
			if (graf2)
				fprintf( out2,"\n  %5.3lf ",f);
			
			//u = sqrt(EPS_MAX);
			u = 1.0;
			for ( m=1; m <=3; m++ )
				X[m] = Abet[m] * u;

			for ( n=1; n <= 21; n++ )
				H_[n-1] = 0.0;
		//__________________________________________________________________
			if (!Nx)//Простой планарный волновод
			{
				hru = 1;
				for ( n=1; n <= 21; n++ )
					H_[n-1] = 0.0;
				m = 1;
				Nsl = Ny + 1;
				
				EPS_r_[1] = eps_sub_r[1];
				EPS_i_[1] =  eps_sub_i[1];

				EPS_r_[Nsl+1] = eps_sub_r[2];
				EPS_i_[Nsl+1] =  eps_sub_i[2];

				for ( n=1; n <= Ny; n++ )
				{
					EPS_r_[n+1] = eps_r[n][m];
					EPS_i_[n+1] = eps_i[n][m];

					H_[n] = B_[n];
				}
				General_COMPL ( Nsl,X, Nkor, fto,bet_r,bet_i);
				
				for ( i=1; i <= n_cor; i++ )
				{
					printf ("\n      %8.5lf   %8.5lf",bet_r[i],bet_i[i]/bet_r[i]);
					fprintf( out,"\n      %8.5lf   %8.5lf",bet_r[i],bet_i[i]/bet_r[i]);
					if (graf2)
						fprintf( out2,"      %8.5lf   %8.5lf",bet_r[i],bet_i[i]/bet_r[i]);
				}
			}//!Nx
			//__________________________________________________________________
			else
				//разбиение по x
				{
					for ( n=1; n <= 21; n++ )
						H_[n-1] = 0.0;
					for ( m=1; m <=Nx+2; m++ )
					{
						//if (m && m !=Nx+1)
						{
							Nsl = Ny + 1;
							EPS_r_[1] = eps_sub_r[1];
							EPS_i_[1] =  eps_sub_i[1];

							EPS_r_[Nsl+1] = eps_sub_r[2];
							EPS_i_[Nsl+1] =  eps_sub_i[2];

							for ( n=1; n <= Ny; n++ )
							{
								EPS_r_[n+1] = eps_r[n][m];//eps[ny][nx]
								EPS_i_[n+1] = eps_i[n][m];

								H_[n] = B_[n];
							}
						}//m && m !=Nx+1

						for ( n=1; n <= Ny+2; n++ )
						{
						//	printf ("\n        %1d    %1d   %8.5lf   %8.5lf   %8.5lf",m,n,EPS_r_[n],EPS_i_[n],H_[n-1]);
							
						}
				
						hru = 1;

						//General_COMPL ( Nsl ,X, Nkor, fto,bet_r,bet_i);
						General_COMPL ( Nsl ,X, 1, fto,bet_r,bet_i);
				

					//	printf ("\n   %1d   %8.5lf   %8.5lf",m,bet_r[1],bet_i[1]/bet_r[1]);//proveril
						/*printf ("\n   %1d   %8.5lf   %8.5lf",m,bet_r[1],bet_i[1]);
						for ( n=1; n <= Ny+2; n++ )
						{
							printf ("\n        %1d   %8.5lf   %8.5lf   %8.5lf",n,EPS_r_[n],EPS_i_[n],H_[n]);
							H_[n] = B_[n];
						}
						*/

						i = 1;
						//for ( i=1; i <= n_cor; i++ )
						{
							mult(bet_r[i],-bet_i[i], bet_r[i],-bet_i[i], &u,&v);

							eps_ef_r[i][m] = u;
							eps_ef_i[i][m] = v;
						}//i

					//	printf ("\n               %1d   bet %8.5lf   %8.5lf     EPS   %8.5lf   %8.5lf",
					//			m,bet_r[i],bet_i[i],eps_ef_r[i][m],eps_ef_i[i][m]);//!!!!!!!!!!!!!!!!!!!!!!!!!!
					}//m

					// ("\n\n ");//!!!!!!!!!!!!!!!!!!!!!!!!!!
					//for ( i=1; i <= n_cor; i++ )
					{
						Nsl = Nx + 1;
						for ( n=1; n <= 21; n++ )
							H_[n-1] = 0.0;

						for ( m=1; m <=Nx+2; m++ )
						{
							EPS_r_[m] = eps_ef_r[i][m];
							EPS_i_[m] = eps_ef_i[i][m];
							
							if (m !=1 && m !=Nx+2)
								H_[m-1] = A_[m-1];

							
						}//m
						 
						//printf ("\n\n ");//!!!!!!!!!!!!!!!!!!!!!!!!!!
						hru = 0;

						//General_COMPL ( Nsl,X, 1, fto,bet_r,bet_i);
						General_COMPL ( Nsl,X, Nkor, fto,bet_r,bet_i);
				
						for ( i=1; i <= n_cor; i++ )
						{
							printf ("\n      %8.5lf   %8.5lf",bet_r[i],bet_i[i]/bet_r[i]);
							fprintf( out,"\n      %8.5lf   %8.5lf",bet_r[i],bet_i[i]/bet_r[i]);
							if (graf2)
								fprintf( out2,"      %8.5lf   %8.5lf",bet_r[i],bet_i[i]/bet_r[i]);
						}
					}//i
				}//разбиение по x
				//__________________________________________________________________
										
		}/* f */ // LOOP over wavelengths f
			
			
		if (graf2)
			fprintf( out2,"\n     fffffffffff");


		if (graf1)
			fprintf( out1,"\n     fffffffffff");
			
	
		
		fflush( out );

		/*tim ( str);*/

		tim ( str);

		fmat_1d( (void*)bet_r,  1, sizeof( double) );
		fmat_1d( (void*)bet_i,  1, sizeof( double) );

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

