
#include "newlib.h"

void BESSEL(double x,double y, int n, double ftol,
			double *J, double *N)
			/* Функции Бесселя
					   1-го рода ( J[0] + i * J[1] ) и
					   2-го рода ( N[0] + i * N[1] )
					   порядка n
			   комплексного аргумента x + i y.
			   Разложение в ряд
			*/
{
	int k;
	double Pi,dvn,dk,nf,rjk,ijk,rz2,iz2,rgk,igk,rzn,izn,rfk,ifk,u,v,to1,to2;

	Pi = 4.0 * atan(1.0);
	ftol *= ftol;
	rz2 = 0.25 * ( x * x - y * y );
	iz2 = 0.5 * x * y;

	rfk = log( 0.25 * (x * x + y * y) ) + 2.0 * 0.577215664902;
	ifk = 2.0 * atan2(y,x);

	dvn = rzn = nf = 1.0;
	izn = 0.0;

	for ( k=1; k <= n; k++ )
	{
		mult( rzn,izn, 0.5 * x,0.5 * y, &rzn,&izn );
		dk = (double)k;
		nf *= dk;
		rfk -= 1.0 / dk;
	}/* k */

	rjk = rzn / nf;
	ijk = izn / nf;

	k = 0;
	J[0] = J[1] = N[0] = N[1] = 0.0;
	do
	{
		J[0] += rjk;
		J[1] += ijk;

		mult( rjk,ijk, rfk,ifk, &u,&v );
		N[0] += u;
		N[1] += v;

		to1 = ( rjk * rjk + ijk * ijk ) / ( J[0] * J[0] + J[1] * J[1]);
		to2 = ( u * u + v * v ) / ( N[0] * N[0] + N[1] * N[1]);

		dk = (double)(k + 1);
		dvn = - 1.0 / ( dk * (dk + n) );

		mult( rjk,ijk, rz2 * dvn,iz2 * dvn, &rjk,&ijk );

		rfk -= 1.0 / dk + 1.0 / (dk + n);

		k++;
	} while( to1 > ftol || to2 > ftol  );

	if ( n )
	{
		dk = ( n == 1) ? 1.0 : nf / ( double)n;
		divide(dk,0.0, rzn,izn, &rgk,&igk);

		for ( k=0; k < n; k++ )
		{
			N[0] -= rgk;
			N[1] -= igk;
			if ( k == n - 1 )
				break;

			dk = (double)(k + 1);
			dk = 1.0 / ( dk * ( - dk + n ) );
			mult( rgk,igk, rz2 * dk,iz2 * dk, &rgk,&igk );
		}/* k */

	} /* n!=0 */

	N[0] /= Pi;
	N[1] /= Pi;
} /* BES */


