
#include "newlib.h"

void Legendr_P( int N, double x, double **Pnm )
{		
	//Вычисляет массив значений присоединненых полиномов Лежандра. 
	//Определение см. Янке. Отличаются от Градштейна знаком при нечетных m. 
	// N - максимальное значение индекса
	// Pnm[n (нижний индекс), m (верхний индекс)]
	// n <= N, m <= N
	
	int n, m;
	double  n1,m1,sq;
	
	for ( n=0; n <= N; n++ )
	for ( m=n+1; m <= N; m++ )
		Pnm[n][m] = 0.0;

	sq = sqrt(1.0 - x * x);
	Pnm[0][0] = 1.0;
	Pnm[1][0] = x;
	for ( n=0; n < N; n++ )
		Pnm[n+1][n+1] = (2.0 * n + 1.0) * sq * Pnm[n][n];
	
	for ( n=1; n < N; n++ )
	{
		n1 = (double)n;
		for ( m=0; m <= n; m++ )
		{
			m1 = (double)m;
			Pnm[n + 1][m] = ( (1.0 + 2.0 * n1) * x * Pnm[n][m] - (n1 + m1) * Pnm[n - 1][m] ) / (1.0 + n1 - m1);
		}//n
	
	
	} /*  m  */

	
}//Legendr_P

void Legendr( int N, double x, double *Pnm)
{		
	//Вычисляет массив значений  полиномов Лежандра. 
	// N - максимальное значение индекса
	// Pn[n]
	// n <= N
	
	int n;
	double  n1;
	
	
	Pnm[0] = 1.0;
	Pnm[1] = x;
	
	for ( n=1; n < N; n++ )
	{
		n1 = (double)n;
		Pnm[n + 1] = ( (1.0 + 2.0 * n1) * x * Pnm[n] - n1 * Pnm[n - 1] ) / (1.0 + n1);
	}//n

}//Legendr



