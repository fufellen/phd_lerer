
#include <stdio.h>
#include <math.h>
#include <stdlib.h>


#define SIGN(x) ( ( (x) < 0.0 ) ? -1 : ( ( (x) == 0.0 ) ? 0 : 1 ) )
extern FILE *out;

double HORKUZ( double a, double b, double fa, double fb,
			double (*funk)(double), double* ff, double tol)
{
	double c,d,e,fc,tol1,xm,p,q,r,s,ffa;

	M2: c = a;
		fc = fa;
		d = b - a;
		e = d;
		ffa = 0.5 * ( fa + fb );
	M3: if ( fabs(fc) >= fabs(fb) )
			goto M4;
		a = b;  fa = fb;
		b = c;  fb = fc;
		c = a;  fc = fa;
	M4: tol1 = fabs(b) * tol;
		xm = 0.5 * (c - b);
		if ( ( fabs(xm) <= tol1 ) || ( fb == 0.0 ) )
			goto M9;
		if ( ( fabs(e) <= tol1) || ( fabs(fa) <= fabs(fb) ) )
			goto M7;
		if ( a != c)
			goto M5;
		s = fb / fa; p = 2.0 * xm * s;
		q = 1.0 - s; goto M6;
	M5: q = fa / fc; r = fb / fc; s = fb / fa;
		p = s * ( 2.0 * xm * q * ( q - r ) - ( b - a ) * ( r - 1.0 ) );
		q = ( q - 1.0 ) * ( r - 1.0 ) * ( s - 1.0 );
	M6: if ( p > 0.0 )
			q = - q;
		p = fabs( p );
		if ( ( 2.0 * p ) >= ( 3.0 * xm * q - fabs( tol1 * q ) ) ||
				p >= fabs( 0.5 * e * q) )
			goto M7;
		e = d; d = p / q; goto M8;
	M7: d = xm;
		e = d;
	M8: a = b; fa = fb;
		if ( fabs(d) > tol1 )
			b = b + d;
		else
			b = b + tol1 * SIGN( xm );
		fb = (*funk)(b);
/*
		printf( "\n    %8.5lf   %11.3e ",b, fb );
		fprintf( out,"\n    %8.5lf   %11.3e ",b, fb );
*/

		if ( fb * SIGN( fc) > 0.0 )
			goto M2;
		else
			goto M3;
	M9: *ff = fb / ffa;
		return ( b );
}

