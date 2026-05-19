
/*  Пpоцедуpа вычисления комплексного опpеделителя матpицы.
 *
 *   Обpащение:  detcompl( n, d , id);
 *               n - поpядок системы;
 *               (D + i * iD) - элементы опpеделителя
								размера[0:n-1,0:n-1];
				 определитель равен exp(d[1][1]) * (d[0][0] + i * id[0][0])
 *
 */

#include "newlib.h"

void detcompl( int n, double ** d , double ** id )
{
	double mn, t, a, it, ia, max, rm, im, u, v;
	int k, j, i;
	mn = 1.0;
	a = 0.0;
	ia = 0.0;
	for ( k=0; k < n; k++ )
	{
		max = 0.0; j = 0;
		for ( i=k; i < n; i++ )
		{
			t = sqrt(d[i][k] * d[i][k] + id[i][k] * id[i][k]);   /*t = d[i][k];*/
			if (  t  >  max )
			/*if ( fabs( t ) > fabs( max ) )*/
			{
				rm = d[i][k];
				im = id[i][k];
				max = t; j = i;
			}
		}  /* end i */
		if ( max )
		{
			if ( j != k )
			{
				mn = - mn;
				for ( i=k; i < n; i++ )
				{
					t = d[j][i];
					d[j][i] = d[k][i];
					d[k][i] = t;
					t = id[j][i];
					id[j][i] = id[k][i];
					id[k][i] = t;
				}  /* end i */
			}  /* end j != k */

			for ( i=k+1; i < n; i++ )
			{
				divide(d[i][k],id[i][k], rm,im, &t,&it);
				/*t = d[i][k] / max;*/
				for ( j=k+1; j < n; j++ )
				{
					/*
					d[i][j] = d[i][j] - t * d[k][j];
					*/
					mult(t, it, d[k][j], id[k][j], &u, &v);
					d[i][j] = d[i][j] - u;
					id[i][j] = id[i][j] - v;
				}
			}  /* end i */

			/*a = a *  d[k][k];*/
			/*mult(a, ia, d[k][k], id[k][k], &u, &v);*/
			u = sqrt( d[k][k] * d[k][k] + id[k][k] * id[k][k] );
			v = atan2(id[k][k], d[k][k]);
			a += log(u);
			ia += v;
		}
		else
		{
			a = 0.0;
			printf( "\n max( d[i][k] ) = 0.0" );
			exit ( 1 );
		}
	}  /* end k */

	d[1][1] = a;
	d[0][0] = mn * cos(ia);
	id[0][0] = mn * sin(ia);
}  /* end of DETcompl */
