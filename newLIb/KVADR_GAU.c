#include "newLib.h"
/*
   Квадратурные узлы и коэффициенты квадратуры наивысшей точности Гаусса для
   интералов на [a,b]. Крылов стр. n <= 10 */

void KVADR_GAU(int m, double a, double b, double *xk,  double *Ak )
{
	int n;
	double A,B;

	A = 0.5 * (b - a);
	B = 0.5 * (b + a);
		 
	UZLI_GAU( m,xk,Ak );

	for ( n=1; n <= m; n++ )
	{
		Ak[n] *= A;
		xk[n] = A * xk[n] + B;		
	} 	
} /*  KVADR_GAU  */


