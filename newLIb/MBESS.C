#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "newlib.h"
#define SIGN(x) ( ( (x) < 0.0 ) ? -1 : ( ( (x) == 0.0 ) ? 0 : 1 ) )
static void besdr( int, double, double, double* );
static double bess( double, double, double, double, int );
static double gamma( double );

static double zn;
static int s1;

void mbess( int i, int ns, int n, double x[], double nu, double **j )
{
	double mu, f, *b;
	int s, s2;

	s2 = (int)floor( nu );
	mu = nu - (double)s2;
	zn = 1.0;

	for ( s1 = 0; s1 <= n; s1++ )
	   if ( x[s1] <= 1.0 )
	   {
			if ( s1 > i )
				zn = - 1.0;
			f = ( zn > 0.0 ) ? exp( - x[s1] ) : 1.0;
			for ( s = 0; s <= ns; s++ )
				j[s][s1] = bess( s + nu, x[s1], 1.0e-7, nu, i ) * f;
        }
		else
        {
			b = amem_1d( ns+3 ,sizeof( double) );
			if ( b == 0 )
				printf("\n No memory for b in MBESS \n");

			zn = SIGN( s1 - i - 0.5 );
			besdr( ns+s2, mu, x[s1], b );
			for ( s = 0; s <= ns; s++)
				j[s][s1] = b[s+s2];

            fmem_1d(b);
        }
  } /* end mbess */

 void besdr( int p, double nu, double z, double *j )
 {
	int m,f;

	if ( z <= 750.0 )
	{
		if ( zn < 0.0 )
			m = ( z > 100.) ? (int)floor( 0.25 * z + 30.0 ):
				(int)floor( 15.0 + z * ( 0.73 - 0.0033 * z ) );
		else
			m = ( z > 100.) ? (int)floor( 1.045 * z + 30.0 ):
				(int)floor( 14.0 + z * ( 1.43 - z * 0.0023 ) );
	}
	else
	{
		printf("\n Z more 750.0\n");
		exit( 0 );
	}
	if ( p + 2 > m )
		m = p + 2;
	{
		double *i, t, a;
		int n, l;

		if ( (i = ( double* )MemCalloc( m+2, sizeof( double ) ) ) == NULL )
		{
			printf("\nNo memory for i in MBESS.\n" );
			exit( 1 );
		}
		i[m+1] = 0.0;
		i[m] = 5.0e-7;
		f = (int)floor( 1.6 + zn / 2.0 );
		l = ( zn < 0.0 ) ? m : (int)floor( 0.5 * m );

		for ( n = m; n >= 1; n--)
			i[n-1] = 2.0 * ( n + nu ) / z * i[n] - i[n+1] * zn;

		if ( nu <= 1e-7 )
		{
			for ( n = 1, a = 0.0; n <= l; n++ )
				a += i[f*n];
			a = 2.0 * a + i[0];
		}
		else
		{
			t = gamma( nu );
			a = nu * i[0] * t;
			for ( n = 1; n <= l; n++ )
			{
				if ( zn < 0.0 )
				{
					a += 2.0 * nu * ( n + nu ) * t * i[n];
					t *= ( n + 2 * nu ) / ( n + 1 );
				}
				else
				{
					t *= ( nu + n - 1 ) / n;
					a += ( nu + 2 * n ) * i[2*n] * t;
				}
			}
			a *= pow( 2.0, nu );
		}
		for ( n = 0; n <= p; n++ )
		{
			double w;

			w = 0.4 - n / 4.0 + floor( ( n + 0.1 ) / 4.0 );
			zn = (f == 2) ? 1.0 : SIGN( w );
			j[n] = i[n] / a * zn;
		}

		MemFree(i);
	}
 } /* end besdr */


 double bess( double n, double x, double eps, double nu, int i )
 {
	double k, term, sum, te;
	int s;

	if ( x <= 1.0e-7 )
	{
		if ( fabs( n - nu ) < 1.0e-8 )
			return( pow( 0.5, nu ) / gamma( 1.0 + nu ) );
         else
			return( 0.0 );
	}
	else
	{
		s = (int)floor( n - nu + 0.5 );
		k = 0.5 * s;
		sum = pow( 0.5, n ) * pow( x,( double )s ) / gamma( 1.0 + n );
		if ( s1 <= i )
		{
			s = (int)floor( k + 0.1 );
			k = 0.5 * s;
			if ( ( k - floor( k + 0.1 ) ) > 0.1 )
				sum = - sum;
		}
		term = sum;
		k = 0.25 * x * x * zn;
		s = 1;
		do
		{
			term *= k / ( s * ( s + n ) );
			sum += term;
			te = fabs( term / sum );
			s++;
		} while ( te > eps );
		return ( sum );
	}
 }  /* end bess */

double gamma( double x)
{
	double t,r,a[21];
	int i;

	a[0] = 1.0;
	a[1] =  0.577215664901532;
	a[2] = -0.655878071520253;
	a[3] = -0.042002635034095;
	a[4] =  0.166538611382291;
	a[5] = -0.042197734555544;
	a[6] = -0.009621971527876;
	a[7] =  0.007218943246663;
	a[8] = -0.001165167591859;
	a[9] = -0.000215241674114;
	a[10] =  0.000128050282388;
	a[11] = -0.000020134854780;
	a[12] = -0.000001250493482;
	a[13] =  0.000001133027231;
	a[14] = -0.205633841e-6;
	a[15] =  0.6116095e-8;
	a[16] =  0.5002007e-8;
	a[17] = -0.1181274e-8;
	a[18] =  0.104342e-9;
	a[19] =  0.7782e-11;
	a[20] = -0.3696e-11;

	t = 1.0 / x;
	if ( x <= 1.0 )
		goto m2;
	do
	{
		t *= x;
		x -= 1.0;
	} while ( x > 1.0 );

	m2:
	r = a[20];
	for ( i = 19; i >= 0; i-- )
		r = a[i] + x * r;

	return( t / r );
}
