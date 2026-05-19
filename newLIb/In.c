
#include "newLib.h"

//вычисление модифицированной функции Бесселя (МФБ) exp(-x)In(x) с помощью интегрального представления. 
//N показыват на сколько порядок квадратуры превосходит порядок МФБ n
// При N=5 абсолютная погрешность порядка 10-9
double In( int n, double x, int N)
{
	int m,M;
	double f,dM,dn,I,Pi;
	
	Pi = 4.0 * atan(1.0);

	M = n + N + (int)x;
	dM = (double)M;
	dn = (double)n;
		
	for ( I = 0.0, m=1; m <= M; m ++ )
	{
		f = Pi / dM * (-0.5 + m);
						
		I += cos(dn * f ) * exp(x * (cos(f) - 1.0));
				
	}//m
		
	I /= dM;
		
	//printf( "\n %3d  %6.3e  %9.8e",n,x,I);
				
	return (I);				

}/* In */

