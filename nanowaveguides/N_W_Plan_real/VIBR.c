
#include "VIBRr.h"

void EPS_MET(int is, double *EPS_r, double *EPS_i);

int main ( void )
{
	int var, Nvar,m;
		//,n,m,M_z1,M_r1,m_z[151],m1,mu,nu,N_g1,mh,nn,mm,N_c;
	double u;
	
	
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
/*//Для ввода параметров материала
	for ( m=200; m <= 2000; m++ )
	{
		double x,n,k;
		fscanf (in,"%lf %lf %lf",&x,&n,&k);
		fprintf ( out,"\n Lam[%1d] =%4.1lf; n[%1d] =%6.3lf;  k[%1d] =%6.3lf; ",m,x,m,n,m,k);
	}
	
	fflush( out );
	exit( 1 );
	*/
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
		
		
		for ( f = AF[1]; f <= AF[3]; f += AF[2] )
		{ 
			k = Pi_2 / f;
			k2 = k * k;

			if (bum1)
			{
				for ( m=1; m <= Nsl+1; m++ )//Ввод диэл. прог. металлических слоев
				{
					if (is[m])
					{
						if ( is[m] == 2 )
						{
							//EPS_r[m] = Ag(f,0);
							//EPS_i[m] = Ag(f,1);
							Ag_(f,&EPS_r[m],&EPS_i[m]);
						}
						else 
							if ( is[m] == 3 )
							{
								//EPS_r[m] = Au(f,0);
								//EPS_i[m] = Au(f,1);
								Au_(f,&EPS_r[m],&EPS_i[m]);
							}
							else
								if ( is[m] == 1 )
								{
									EPS_r[m] = Cu(f,0);
									EPS_i[m] = Cu(f,1);
								}

						EPS_r_[m] = EPS_r[m];
						EPS_i_[m] = EPS_i[m];
					}//is[n]
					
				}//m
			}//bum1
			

			EPS_min = EPS_r[1];
			if (Nsl > 1)
				if (!is[Nsl+1] && EPS_min > EPS_r[Nsl+1])
					EPS_min = EPS_r[Nsl+1];
			
			//Нахождение EPS_MAX
				
			EPS_MAX = tgd = 0.0;
			
			for ( m=1; m <= Nsl+1; m++ )
			{ 	
				u = fabs(EPS_i[m] / EPS_r[m]);
				if (tgd < u)
					tgd = u;
		
				if (is[m] ||  EPS_r[m] < 0.0)
					continue;
				
				u = fabs(EPS_r[m]);
				if (EPS_MAX < u)
					EPS_MAX = u;
			} 

			//Нахождение EPS_MAX END
			
			printf ("\n\n  lamda =  %5.3lf ",f);
			fprintf( out,"\n\n  lamda =  %5.3lf ",f);
			if (graf2)
				fprintf( out2,"\n  %5.3lf ",f);
			
			//u = sqrt(EPS_MAX);
			u = 1.0;
			for ( m=1; m <=3; m++ )
				X[m] = Abet[m] * u;
			
			for ( m=4; m <=5; m++ )
				X[m] = Abet[m] * u * tgd;
			
		
			if (TV == 3)
				F_1();
			else
				General_COMPL ( X, Nkor, fto);
							
		}/* f */ // LOOP over wavelengths f
			
			
			if (graf2)
				fprintf( out2,"\n     fffffffffff");


			if (graf1)
				fprintf( out1,"\n     fffffffffff");
			
	
		
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