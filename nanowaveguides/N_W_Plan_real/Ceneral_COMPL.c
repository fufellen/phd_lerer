
#include "VIBRr.h"
double GENERAL(double COM,double *X,double func(double f),double *ff );
double F_df( double x);

//void mult_a( double *x,double *y,double *z );
//void div_a( double *x,double *y,double *z );
void subt_a( double *x,double *y,double *z );
void addit_a( double *x,double *y,double *z );
double mod( double *x,double *y);
void NEW_C( double *x, double st, void func(double *z,double *f), double* ff, double tol);
int sig(double x);

int FIN;
double zn,FIN1;

void General_COMPL ( double *X, int N_kor,double fto )
{
	//int N_c,max=1000;
	int i,n_cor,m,M,n;
	double COM,F,ff;
		
	/********************
	// Комплексные корни 29.11.12. Тест на sin(a*Pi*z)=0, a=a1+i*a3
	Вначале ищется решение уравнения sin(a1z)=0. Решение с помощью нашей процедуры поиска корня.
	Затем последовательно методом Ньютона ищутся решения уравнений sin((a1+i*a2) * z)=0.
	Величина a2 меняется от 0 с шагом st2 от до a3. 
	Полученное приближенное решение на к-м шаге - начальное значения для шага к+1.
	fto - относительная погрешность
	*******************/

	zn = (X[3] > X[1]) ? 1.0 : -1.0;
	
	FIN1 = sqrt(EPS_min) + fto;
	if (X[3] > FIN1)  
		FIN1 = X[3];
	//Вначале полагаем, что потерь нет и ищем действительное решение.
	//Для волновода без потерь это строгое окончательное решение
	//Для волновода с потерями для контроля комплексного решения
	
	M = (Pot) ? 3 : 1;

	for ( m=1; m <= M; m++ )
	{
		ff = (m == 1) ? 0.0 : (m == 2) ? 1.0 : -1.0;
		
		for ( n=1; n <= Nsl+1; n++ )
		{
			EPS_r[n] =  EPS_r_[n] + ff * EPS_i_[n];
			EPS_i[n] = 0.0;
		}

		COM = X[1];

		for ( i=1; i <= Nkor; i++ )
		{	
			ik = 0;
		
			DK[m][i] = F = GENERAL(COM,X,F_df, &ff);

			if ( FIN )
				break;
		
			if (!Pot)
			{
				printf ("\n      %8.5lf   %5.3e ",F,ff);
				fprintf( out,"\n      %8.5lf   %5.3e ",F,ff);

				if (graf2)
					fprintf( out2,"      %8.5lf ",F);
			}
		
			COM = F + zn * 10.0 * fto;
					
		}//i
		if (m == 1)
			n_cor = i - 1;
	}//m
		
	//END  потерь нет .

	
	if (!Pot)
		return;

	printf ("\n  ");
	fprintf( out,"\n  ");
	
	// решение для волновода с потерями
	for ( i=1; i <= n_cor; i++ )
	{	
		double z[3];
		
		ik = 0;
		z[2] = DK[1][i] - DK[2][i];
		z[0] = DK[2][i] + DK[3][i] - 2.0 * DK[1][i];
		z[1] = DK[1][i] - z[0] * 0.5;
		
		//ff = z[1] / sqrt(EPS_min);
		ff = z[1] / sqrt(EPS_MAX);
		
		printf ("\n      %8.5lf   %8.5lf   %8.5lf",z[1],z[2]/z[1],ff);
		fprintf( out,"\n      %8.5lf   %8.5lf   %8.5lf",z[1],z[2]/z[1],ff);
		if (graf2)
			fprintf( out2,"      %8.5lf   %8.5lf   %8.5lf",z[1],z[2]/z[1],ff);
					
	}//END  потерь нет .

}   //General_COMPL

//_____________________________________________________________

double GENERAL(double COM,double *X,double func(double f),double *f )
{
	double fe,a,b,fa,fb,ff;

	FIN = 0;
	a = COM;
		
	fa = F_df(a);

	met: b = a + zn * X[2];
	
	if ( (zn > 0.0 && b > FIN1 * 1.000001) ||  (zn < 0.0 && b < FIN1 * 1.000001)   )
	{
		b = FIN1 * 1.000001;
		FIN = 1;
	}
	fb = F_df(b);
	
	
	if ( sig(fa) == sig(fb) )
	{
		a = b;
		fa = fb;
		if (!FIN)
			goto met;
		else
			return (0.0);
	}
	else
		fe = HORKUZ( a, b, fa, fb, F_df, &ff, fto );
		
	*f = F_df(fe);
	
	return (fe);

}   /*  GENERAL  */

//_____________________________________________________________


int sig(double x)
{
	int y;
	y = ( x > 0.0 ) ? 1 : ( x < 0.0 ) ? -1 : 0;
	return (y);
}

double F_df( double x) 
{
	double y,f[3],z[3];		
	if( kbhit() )
	   if( getch() == ESCape )
		 {
		 fflush( out );
		 exit( 1 );
		 }

	z[1] = x;
	z[2] = 0.0;

	F_Compl(z,f);//!!!!!!!!!!!!!!!!!!!!!!!!!!!Изменить

	y = (fabs(f[1]) > 1.0e-31) ? f[1] : f[2];

//	printf ("\n     %5.3e  %5.3e    %5.3e  %5.3e",z[1],z[2],f[1],f[2]);
//	fprintf( out,"\n     %5.3e  %5.3e    %5.3e  %5.3e",z[1],z[2],f[1],f[2]);


	return (y);	
		
}/* F */


void mult_a( double *x,double *y,double *z )//компл умножение
{
	z[1] = x[1] * y[1] - x[2] * y[2];
	z[2] = x[1] * y[2] + y[1] * x[2];
}

void div_a( double *x,double *y,double *z )//компл деление
{
	double zz;
	zz = 1.0 / ( y[1] * y[1] + y[2] * y[2] );
	z[1] = (x[1] * y[1] + x[2] * y[2]) * zz;
	z[2] = (-x[1] * y[2] + y[1] * x[2]) * zz;
}

void subt_a( double *x,double *y,double *z )//компл вычитание
{
	z[1] = x[1] - y[1];
	z[2] = x[2] - y[2];
}

void addit_a( double *x,double *y,double *z )//компл сложение
{
	z[1] = x[1] + y[1];
	z[2] = x[2] + y[2];
}

double mod( double *x,double *y)//модуль разности двух компл чисел
{
	double xx[3];
	xx[1] = x[1] - y[1];
	xx[2] = x[2] - y[2];
	return (sqrt(xx[1] * xx[1] + xx[2] * xx[2]));
}