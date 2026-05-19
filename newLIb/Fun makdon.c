#include "newLib.h"



double K0( double x )
{
	/* Функция Макдональда */
	double z,I,t;
	if ( x <= 2.0 )
	{
		z = x / 3.75;
		z *= z;
		I = 1.0 + z * ( 3.51562 + z * ( 3.08994 + z * ( 1.20674 +
			z * ( 0.26597 + z * ( 0.03607 + z * 0.00458 )))));
		t = x / 2.0;
		z = t * t;
		z = ( - log( t ) * I - 0.57721566 + z * ( 0.42278420 + z * ( 0.23069756
			+ z * ( 0.03488590 + z * ( 0.00262698 + z * ( 0.00010750
			+ z * 0.00000740 ))))));
	}
	else
	{
		z = 2.0 / x;
		z = ( 1.25331414 + z * ( - 0.07832358 + z * ( 0.02189568
			+ z * ( - 0.01062446 + z * ( 0.00587872 + z * ( - 0.00251540
			+ z * 0.00053208 ))))) )* exp(-x);
		z /= sqrt(x);
	}
	return (z);
} /* end K0 */


double K1(double x)
{
	double t, t2, z;
	t = x * 0.5;
	t2 = t * t;
	if (x < 2.0)
	{
		z = x * log(t) * I1(x) + 1.0 + t2 * (0.15443144 + t2 * (-0.67278579
					+ t2 * (-0.18156897 + t2 * (-0.01919402 + t2 *
					(-0.00110404 - t2 * 0.00004686)))));
							/* precision = 8e-9 */
		z /= x;
	}
	else
	{
		t = 2.0 / x;
		z = 1.25331414 + t * (0.23498619 + t * (-0.03655620 + t * (0.01504268
			+ t * (-0.00780353 + t * (0.00325614 - t * 0.00068245)))));
		   /* precision = 2.2e-7*/
		z *= exp(-x) / sqrt(x);
	}
	return(z);
} /* end K1 */

double Kn( int N,double x )/* Функция Макдональда */
{
	int n;
	double k,km,kp;
	if ( !N )
		kp = K0(x);
	else
		if ( N == 1 )
			kp = K1(x);
		else
		{
			k = K1(x);
			km = K0(x);
			for ( n=1; n < N; n++ )
			{
				kp = km + 2.0 * n * k / x;
				km = k;
				k = kp;
			}
		}
	return(kp);
}/* end Kn */


double Knp( int N,double x )/* Производная от Функции Макдональда */
{
	int n;
	double k,km,kpp,kp;
	if ( !N )
		kpp = - K1(x);
	else
	{
		k = K1(x);
		km = K0(x);
		for ( n=1; n <= N; n++ )
		{
			kpp = - km - k * n / x;
			kp = - km - 2.0 * kpp;
			km = k;
			k = kp;
		}
	}
	return(kpp);
} /* end Knp */


