
/*  Пpоцедуpа вычисления знака опpеделителя матpицы.
 *
 *   Обpащение:  detsig( n, D );
 *               n - поpядок системы;
 *               D[0:n-1,0:n-1] - элементы опpеделителя;
 *
 */

#include "newlib.h"

#define SIGN(x) ( ( (x) < 0.0 ) ? -1 : ( ( (x) == 0.0 ) ? 0 : 1 ) )
//#include <matrix.h>
//#include <lib.h>

double detsig( int n, double **d )
{
	double t, a, max, nd;
	int k, j, i;
	a = 1.0; nd = 0.0;
	for ( k=0; k < n; k++ )
	{
		max = 0.0; j = 0;
		for ( i=k; i < n; i++ )
		{
			t = d[i][k];
			if ( fabs( t ) > fabs( max ) )
			{
				max = t; j = i;
			}
		}  /* end i */
		if ( max )
		{
			if ( j != k )
			{
				a = - a;
				for ( i=k; i < n; i++ )
				{
					t = d[j][i];
					d[j][i] = d[k][i];
					d[k][i] = t;
				}  /* end i */
			}  /* end j != k */

			for ( i=k+1; i < n; i++ )
			{
				t = d[i][k] / max;
				for ( j=k+1; j < n; j++ )
					d[i][j] = d[i][j] - t * d[k][j];
			}  /* end i */
			a = a * SIGN( d[k][k] );
			nd += 0.434294482 * log( fabs( d[k][k] ) );
		}
		else
		{
			a = 0.0;
			printf( "\n max( d[i][k] ) = 0.0" );
			exit ( 1 );
		}
	}  /* end k */
	d[0][0] = nd;

	return ( a );

}  /* end of DETSIG */
