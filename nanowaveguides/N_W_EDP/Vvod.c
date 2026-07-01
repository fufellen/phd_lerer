
#include "VIBRr.h"

void VVOD ( void )
{
	int n,m;
	double n_,k_;
	
	
	// Нановолноводы. Метод ЭДП. 16.10.13
	
	
	/*
	Параметры цикла по длине волны AF
	Длина волны и размеры а нанометрах
	*/

	fSkipComment( in );
	fscanf (in,"%lf %lf %lf",&AF[1],&AF[2],&AF[3]);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	
	/*
	E = 1 - E-волна. В плазмонных волноводах может распространяться только эта волна
	E = 0 - H-волна.
	
		Если is[n] = 0, то считывается из файла с ч.д
		{
			n_[n], k_[n] - действительная и мнимая части показателя преломления.
		}
		Если is[n] = 1, то n_[n], k_[n] определяются в базе данных для металлов.
		is[n]= 1 медь
		is[n]= 2 - серебро
		is[n]= 3 - золото

		
	 */
	fSkipComment( in );
	fscanf (in,"%d %d %d",&E,&Nx_,&Ny);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	if (Nx_ > 1)
	{
		Nx = Nx_ - 2;
		nsub = 2;
		nc = 1;
	}
	else
	{
		Nx = 0;
		nsub = 2;
		nc = 1;
	}
	
	Pot = 0;
	for (  bum1 = 0,n=1; n <= Ny; n++ )
	{
		fSkipComment( in );
		for ( m=1; m <= Nx_; m++ )
		{ 
			fscanf (in,"%d",&is[n][m]);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
			bum1 += is[n][m];
	
			if (!is[n][m])
			{
				fscanf (in,"%lf %lf ",&n_,&k_);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				
				//eps_r_[n][m] = eps_r[n][m] = n_[n] * n_[n] - k_[n] * k_[n];
				//eps_i_[n][m] =  eps_i[n][m] =  - 2.0 * n_[n] * k_[n];
				eps_r_[n][m] = eps_r[n][m] = n_ * n_ - k_ * k_;
				eps_i_[n][m] =  eps_i[n][m] =  - 2.0 * n_ * k_;
					
			}//!is[n]
					
			if (is[n][m] || fabs(eps_i[n][m]) > 1.0e-6)
				Pot = 1;

		} //m 
	}//n

	
	for (  bum2 = 0,n=1; n <= 2; n++ )
	{ 
		fSkipComment( in );
		fscanf (in,"%d",&is_sub[n]);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				
		bum2 += is_sub[n];
	
		if (!is_sub[n])
		{
			fscanf (in,"%lf %lf ",&n_,&k_);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				
			eps_sub_r[n] = n_ * n_ - k_ * k_;
			eps_sub_i[n] =  - 2.0 * n_ * k_;
					
		}//!is[n]
					
		if (is_sub[n]|| fabs(eps_sub_i[n]) > 1.0e-6)
			Pot = 1;

	} //n 
		
	
		for ( m=1; m <= Nx; m++ )
		{
			fSkipComment( in );
			fscanf (in,"%lf",&A_[m]);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		}
	
		fSkipComment( in );
		for ( m=1; m <= Ny; m++ )
			fscanf (in,"%lf",&B_[m]);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	/*
	Nkor - число корней.
	fto - относительная точность поиска корня
	*/
	fSkipComment( in );
	fscanf (in,"%d %lf",&Nkor,&fto);
	
	/*
	bet - комплексный коэффициент замедления
	Решения для действительной части корня ищется от Abet[1] * sqrt(EPS_MAX)   до Abet[3] * sqrt(EPS_MAX) !!! Abet[3]< Abet[1]
	с шагом Abet[2] * sqrt(EPS_MAX).
	
	EPS_MAX - наибольшее из всех епселонов по модулю
	tgd - наибольшее из всех тангенс делта по модулю
	Следует учесть, что для обычного диэлектрического волновода, в том числе с потерями, Abet[1] < 1.
	Для плазмонного волновода может быть и Abet[1] > 1.

	
	*/
	fSkipComment( in );
	fscanf (in,"%lf %lf %lf",&Abet[1],&Abet[2],&Abet[3]);

//Конец ввода числовых данных

  	
//____________________________________________________________________________________________________________	
	if (E)
	{
		printf ( "\n  E-wave. ");
		fprintf ( out,"\n  E-wave.");
	}
	else
		{
			printf ( "\n  H-wave.");
			fprintf ( out,"\n  H-wave. ");
		}

	if (Razb)
	{
		printf ( "\n  Razbienie po x ");
		fprintf ( out,"\n  Razbienie po x ");
	}
	else
		{
			printf ( "\n  Razbienie po y ");
			fprintf ( out,"\n  Razbienie po y ");
		}

	printf ( "\n\n  Waveguide ");
	fprintf ( out,"\n\n  Waveguide ");
	
	for ( n=1; n <= Ny; n++ )
	{
		printf ( "\n  ");
		fprintf ( out,"\n  ");
		for ( m=1; m <= Nx_; m++ )
		{ 
			if (is[n][m])
			{
				printf ( "\n  n=%1d   m=%1d metall is[n,m]=%1d",n,m,is[n][m]);
				fprintf ( out,"\n  n=%1d   m=%1d metall is[n,m]=%1d",n,m,is[n][m]);
			}
			else
			{
				printf ( "\n  eps_r[%1d,%1d]=%8.5lf  eps_i[%1d,%1d] =%8.5lf",
					n,m,eps_r[n][m],n,m,eps_i[n][m]);
				fprintf ( out,"\n  eps_r[%1d,%1d]=%8.5lf  eps_i[%1d,%1d] =%8.5lf",
					n,m,eps_r[n][m],n,m,eps_i[n][m]);
			}

		} //m 
	}//n
	printf ( "\n  ");
	fprintf ( out,"\n  ");
	for ( m=1; m <= Nx; m++ )
	{
		printf ( "\n  m=%1d   A_[m]=%8.5lf",m,A_[m]);
		fprintf ( out,"\n  m=%1d   A_[m]=%8.5lf",m,A_[m]);
	}
	printf ( "\n  ");
	fprintf ( out,"\n  ");			
	for ( n=1; n <= Ny; n++ )
	{
		printf ( "\n  n=%1d   B_[n]=%8.5lf",n,B_[n]);
		fprintf ( out,"\n  n=%1d   B_[n]=%8.5lf",n,B_[n]);
	}

	printf ( "\n\n  Substrate");
	fprintf ( out,"\n\n  Substrate");

	for ( n=1; n <= nsub; n++ )
	{ 
		if (is_sub[n])
		{
			printf ( "\n  n=%1d   metall is[n]=%1d",n,is_sub[n]);
			fprintf ( out,"\n  n=%1d   metall is[n]=%1d",n,is_sub[n]);
		}
		else
		{
			printf ( "\n  eps_sub_r[%1d]=%8.5lf  eps_sub_i[%1d] =%8.5lf",n,eps_sub_r[n],n,eps_sub_i[n]);
			fprintf ( out,"\n  eps_sub_r[%1d]=%8.5lf  eps_sub_i[%1d] =%8.5lf",n,eps_sub_r[n],n,eps_sub_i[n]);
		}

	} // 
		
	
		printf ( "\n\n  Nkor=%2d   fto=%5.3e",Nkor,fto);
		fprintf ( out,"\n\n  Nkor=%2d   fto=%5.3e",Nkor,fto);
		for ( n=1; n <= 3; n++ )
		{
			printf ( "     X[%1d]=%4.2e ",n,Abet[n]);
			fprintf ( out,"     X[%1d]=%4.2e ",n,Abet[n]);
		}
		
}   /*  end VVOD  */


	  
