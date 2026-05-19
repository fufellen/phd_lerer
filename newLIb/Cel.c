/*
 *   Процедура решения системы линейных уравнений с действительными
 *      коэффициентами
 *   Обращение: CEL( n, nm, a );
 *        n - ранг матрицы ( число строк)
 *        nm = n + m, m - число столбцов правой части
 *        a[1:n,1:nm] - действительная матрица размера n * nm
 */
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

void CEL( int n, int nm, double **a)
{
	int k, l, i, i1, j;
	double d, br, za, zc, zb;

	for ( k=1; k <= n; k++ )
	{
		d = 0.0;
		l = 1;
		for ( j=k; j <= n; j++ )
		{
			br = fabs( a[j][k] );
			if ( br > d)
			{
				d = br;
				l = j;
			} /* end if */
		 }  /* end j */
		 if ( d < 1.0e-19 )
		 {
			printf("\n CEL: STOP!  max( fabs( a[i,j] ) ) = 0.0 \n");
			exit(1);
		 }
		 za = a[l][k];
		 a[l][k] = a[k][k];
		 for ( j=k+1; j <= nm; j++ )
		 {
			zc = a[l][j];
			if ( d < fabs( zc ) * 1.0e-19 )
			{
				printf("\n CEL : STOP! \n");
				exit( 1 );
			}
			a[l][j] = a[k][j];
			zb = zc / za;
			a[k][j] = zb;
			for ( i=k+1; i <= n; i++ )
			{
			   a[i][j] -= zb * a[i][k];
			}
		 }
	}  /*  end k */

	for ( l=n+1; l <= nm; l++ )
	{
		for ( i1=1; i1 <= n-1; i1++ )
		{
			i = n - i1;
			za = a[i][l];
			for ( j=i+1; j <= n; j++ )
			   za -= a[j][l] * a[i][j];
			a[i][l] = za;
		}  /*  end i1 */
	}  /*  end l */
}
