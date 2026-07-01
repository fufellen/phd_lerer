
#include "VIBRr.h"

void VVOD ( void )
{
	int n;
	
	
	// Планарные нановолноводы. 06.01.13
	
	
	/*
	Параметры цикла по длине волны AF
	Длина волны и размеры а нанометрах
	*/

	fSkipComment( in );
	fscanf (in,"%lf %lf %lf",&AF[1],&AF[2],&AF[3]);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	
	/*
	E = 1 - E-волна. В плазмонных волноводах может распространяться только эта волна
	E = 0 - H-волна.

	TV  управляющая переменная:
		а) TV = 1 или TV = 2, то это волна в симметричном волноводе. Слой с номером 2 имеет толщину  2 *H_[2].
				TV = 1 в плоскости симметрии магнитная стенка, TV = 2 в плоскости симметрии электрическая  стенка
		
		б) TV = 5, TV = 3,то это волна на границе двух полубесконечных слоев. Nsl = 1. 
			Корень находится аналитически при TV = 3;
			Корень находится численно (для тестирования процедуры поиска комплексного корня) при TV = 5.
		
		в) TV = 4  - волновод из Nsl + 1 слоев. Нумерация диэлектриков сверху вниз. Общее число слоев  Nsl + 1.
		pervii polubeskonechnii, слой с номером Nsl + 1 - polubeskonechnaia podlojka.
		Если TV = < 4, т.е. варианты а),б) Nsl = 1.		
		
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
	fscanf (in,"%d %d %d",&E,&TV,&Nsl);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	Pot = 0;
	for ( bum1 = 0, n=1; n <= Nsl+1; n++ )
	{ // loop over 6 materials
		//if (n == 1)
		//	is[n] = 0;
		//else
		{
			fSkipComment( in );
			fscanf (in,"%d",&is[n]);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		}
		bum1 += is[n];
	
		if (!is[n])
		{
			if (n ==1 || (n == Nsl+1 && TV > 2))
			{
				fSkipComment( in );
				fscanf (in,"%lf %lf ",&n_[n],&k_[n]);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				H_[n-1] = 0.0;
			}//(n ==1 || n == Nsl+1
			else
			{				
				fscanf (in,"%lf %lf %lf ",&H_[n-1],&n_[n],&k_[n]);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			}//(n !=1 && n != Nsl+1
				
					
			EPS_r_[n] = EPS_r[n] = n_[n] * n_[n] - k_[n] * k_[n];
			EPS_i_[n] =  EPS_i[n] =  - 2.0 * n_[n] * k_[n];
					
		}//!is[n]
		else
			if (!(n ==1 || (n == Nsl+1 && TV > 2)))
			{
				fscanf (in,"%lf ",&H_[n-1]);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				
			}//is[n]
		
		if (is[n] || fabs(EPS_i[n]) > 1.0e-6)
			Pot = 1;

	} //n // loop over 6 materials

		
		
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
	
	if (TV == 1)
	{
		printf ( "  Magnitnaya stenka ");
		fprintf ( out,"  Magnitnaya stenka");
	}
	else
		if (TV == 2)
		{
			printf ( "  Electricheskaya stenka");
			fprintf ( out,"  Electricheskaya stenka ");
		}

	for ( n=1; n <= Nsl+1; n++ )
	{
		if (n > 1 && (n < Nsl+1 || TV < 3))
		{
			if (is[n])
			{
				printf ( "\n  n=%1d metall is[n]=%1d  H[%1d] =%6.3lf",n,is[n],n,H_[n-1]);
				fprintf ( out,"\n  n=%1d metall is[n]=%1d  H[%1d] =%6.3lf",n,is[n],n,H_[n-1]);
			}
			else
			{
				printf ( "\n  n[%1d]=%8.5lf  k[%1d] =%8.5lf  H[%1d] =%6.3lf",n,n_[n],n,k_[n],n-1,H_[n-1]);
				fprintf ( out,"\n  n[%1d]=%8.5lf  k[%1d] =%8.5lf  H[%1d] =%6.3lf",n,n_[n],n,k_[n],n-1,H_[n-1]);
			}
		
		}
		else
		{
			if (is[n])
			{
				printf ( "\n n=%1d metall is[n]=%1d",n,is[n]);
				fprintf ( out,"\n n=%1d metall is[n]=%1d",n,is[n]);
			}
			else
			{
				printf ( "\n  n[%1d]=%8.5lf  k[%1d] =%8.5lf",n,n_[n],n,k_[n]);
				fprintf ( out,"\n  n[%1d]=%8.5lf  k[%1d] =%8.5lf",n,n_[n],n,k_[n]);
			}
						
		}
	}//n

		printf ( "\n\n  Nkor=%2d   fto=%5.3e\n",Nkor,fto);
		fprintf ( out,"\n\n  Nkor=%2d   fto=%5.3e\n",Nkor,fto);
		for ( n=1; n <= 3; n++ )
		{
			printf ( "X[%1d]=%4.2e ",n,Abet[n]);
			fprintf ( out,"X[%1d]=%4.2e ",n,Abet[n]);
		}
		
}   /*  end VVOD  */


	  
