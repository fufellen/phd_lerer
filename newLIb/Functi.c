
//#include "cylindr.h"

  /* Bessel's functions */

/* modified Bessel function of zero order */    double I0(double x);
/* modified Bessel function of first order */   double I1(double x);
/* Mc'Donald function of zero order */          double KKK0(double x);
/* Mc'Donald function of 1-th order */          double KKK1(double x);
//Процедуа вычисления производной модифицированных
// функций Бесселя Im(r1). Метод разложения в ряд
//Процедуа вычисления модифицированных функций Бесселя 
//Im(r1). Метод разложения в ряд
#include "newlib.h"


double I0(double x)
{
	double t, KK;
	if ( x < 3.75 )
	{
		t = (x/3.75)*(x/3.75);
		KK = 1.0 + t*( 3.5156229 + t*(3.0899424 + t*(1.2067492 +
				 t*(0.2659732 + t*(0.0350768 + t*0.0045813 )))));
						  /* precision 1.6e-8 */
	}
	else
	{
		t = 3.75/x;
		KK = 0.39894228 + t*(0.01328592 + t*( 0.00225319 + t*(-0.00157565 + t*(0.00916281 +
				  t*(-0.02057706 + t*(0.02635537 +
				  t*(-0.01647633 + t*0.00392377)))))));
							/* precision 1.9e-7 */
		KK= KK/(sqrt(x))*exp(x);
	}
	return(KK);
}

/* modified Bessel function of first order */
double I1(double x)
{
	double t, KK;
	if ( x <   3.75 )
	{
		t = (x / 3.75) * ( x / 3.75);
		KK = 0.5 + t*(0.87890594 + t*(0.51498869 + t*(0.15084934 +
			   t*(0.02658733 + t*(0.00301532 + t*0.00032411)))));
							   /* presicion 8e-9 */
		KK = KK*x;
	}
	else
	{
		t = 3.75/x;
		KK = 0.39894228 + t*(-0.03988024 + t*(-0.00362018 + t*(0.00163801 +
						  t*(-0.01031555 + t*(0.02282967 + t*(-0.02895312 +
						  t*(0.01787654 - t*0.00420059)))))));
										 /* presicion 2.2e-7 */
		KK = KK*exp(x)/sqrt(x);
	}
	return(KK);
}

double KKK0(double x)
{
	double t,KK;
	if (x< 2.0)
	{
		t = (x/ 2.0)* (x/ 2.0);
		KK= - log(x/ 2.0)* I0(x) - 0.57721566 + t*(0.42278420 +
							   t*(0.23069756 + t*(0.03488590 +
							   t*(0.00262698 + t*(0.00010750 +
							   t*0.00000740)))));
							   /* precision 1e-8 */
	}
	else
	if (x>= 2.0)
	{
		t =  2.0/x;
		KK = 1.25331414 + t*(-0.07832358 + t*(0.02189568 + t*(-0.01062446 +
					   t*(0.00587872 + t*(-0.00251540 + t*0.00053208)))));
						  /* precision 1.9e-7 */
		KK = KK /sqrt(x)*exp(-x);
	}
	return(KK);

}

double KKK1(double x)
{
	double t, k,  KK;
	  t = x/ 2.0;
	  k = t * t;
	
	if (x< 2.0)
	{
		  KK = x*log(t)*I1(x) + 1.0 + k*(0.15443144 + k*(-0.67278579 + k*(-0.18156897 +
			k*(-0.01919402 + k*(-0.00110404 - k*0.00004686)))));
							  /* precision = 8e-9 */
		  KK = KK/x;
	}
	else
	if (x>= 2.0)
	{
	 t =  2.0/x;
	 KK = 1.25331414 + t*(0.23498619 + t*(-0.03655620 + t*(0.01504268 +
		 t*(-0.00780353 + t*(0.00325614 - t*0.00068245)))));
		/* precision = 2.2e-7*/
	 KK = exp(-x)*KK/sqrt(x);
	 }
  return(KK);
}



//Процедуа вычисления модифицированных функций Бесселя Im(r1). Метод разложения в ряд

double Im( int m, double z )
{
	double kd,t,t2,jk,I;
	int k;

	t = 0.5 * z;
	t2 = t * t;

	for ( jk = 1.0, k=0; k < m; k++ )
		jk *= t / (1.0 + k);
	
	I = jk;
	k = 0;
	do
	{
		k++;
		kd = (double)k;
		jk *= t2 / (kd * (kd + m));
		I += jk;
	} while( jk / I >  1.0e-7 );
	
	
	return(I);
	
}/* Im */


//Процедуа вычисления производной модифицированных функций Бесселя Im(r1). Метод разложения в ряд

double Imp( int m, double z )
{
	double kd,t,t2,jk,jk1,I;
	int k;
	
	if ( z < 1.0e-12)
		return(0.0);

	t = 0.5 * z;
	t2 = t * t;

	for ( jk = 1.0, k=0; k < m; k++ )
		jk *= t / (1.0 + k);
	
	I = jk * m;
	k = 0;
	do
	{
		k++;
		kd = (double)k;
		jk *= t2 / (kd * (kd + m));
		jk1 = jk * (2.0 * kd + m);
		I += jk1;
		if (k == 1)
			jk1 = 1.0e10;
	} while( jk1 / I >  1.0e-7 );
	
	
	return(I / z);
	
}/* Imp */

