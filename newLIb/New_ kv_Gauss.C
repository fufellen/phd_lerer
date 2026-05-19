#include "newLib.h"
/*
   Квадратурные узлы и коэффициенты квадратуры наивысшей точности Гаусса для
   интералов на [a,b]. Число квадратурных узлов M  - произвольное */


void NEW_KV_GAUS(int M, double a, double b, double *xk,double *Ak )
{
	double an,bn,d,xk2[11],Ak2[11];
	int N,p,n,M1;
	
	if (M < 11)
		KVADR_GAU(M, a, b, xk,  Ak );
	else
	{
		N = (int)(M / 10);
		d = (b - a ) * 10.0 / M;
	
		//x = ( double**)amat_2d( 1,N, 1,10, sizeof( double) );
	//	A = ( double**)amat_2d( 1,N, 1,10, sizeof( double) );

		for ( n=0; n <= N; n++ )
		{
			an = a + d * n;
			if ((n != N))
			{
				M1 = 10;
				bn = an + d;
			}
			else
			{
				M1 = M - 10 * N;
				bn = b;
			}
		
			KVADR_GAU(M1, an, bn, xk2,  Ak2 );
			
			for ( p=1; p <= M1; p++ )
			{
				xk[10*n+p] = xk2[p];
				Ak[10*n+p] = Ak2[p];
			}
		}
		
		
	//	fmat_2d( (void**)x, 1,N, 1, sizeof( double) );
	//	fmat_2d( (void**)A, 1,N, 1, sizeof( double) );
	}

		
}   /* NEW_KV_GAUS */


