
/*  пpоцедуpа pешения системы линейных алгебраических уpавнений (Xj * Amj)=Bm
 *    с комплексными коэффициентами
 *      ( с двойной точностью )
 *   СЛАУ приводится к виду: C Xm+(Xj * Amj)=Bm, 
 *		констаета C задается в списке параметров процедуры.
 *
 *   Обpащение:  ITRSLAU( M, I, ftol,RA, IA, ReC, ImC );
 *               M - порядок СЛАУ;
 *               I - максимальное количество иттераций;
				 ftol - относительная точность решения СЛАУ;
 *               RA - действительная часть коэффициентов;
 *               IA - мнимая часть коэффициентов;
 *               RA[1:M,1+M] - действительная часть коэффициентов в правой части;
 *               IA[1:m,1+M] - мнимая часть коэффициентов в правой части;
 *               ReC - действительная часть констаеты C;
 *               ImC - мнимая часть констаеты C;
 *				 V[1:M]	вес
 *					Решение в  матричных элементах
 *               RA[1:M,1+M] - действительная часть решения;
 *               IA[1:M,1+M] - мнимая часть решения;
 */


#include "newlib.h"

void ITERSLAU(int M, int It, double ftol, double **RA, double **IA,double ReC, double ImC, double *V )
{
	int m,n,i;
	double RS,IS,S1,S,BI,BR,Re,Im,*RX,*IX,Rx,Ix,D;
	
	RX = amat_1d( 1,M, sizeof( double) );
	IX = amat_1d( 1,M, sizeof( double) );
		
	for( m = 1; m <= M; m++)
	{
		for( n = 1; n <= M+1; n++)
		{
			Re = RA[m][n];
			Im = IA[m][n];
			D = 1.0 / ( ReC * ReC + ImC * ImC );
			RA[m][n] = ( Re * ReC + Im * ImC ) * D;
			IA[m][n] = ( Im * ReC - Re * ImC ) * D;
		}
		RA[m][m] -=  1.0;
	}

	//for( m = 1; m <= M; m++)
	//	RA[m][m] -=  1.0;
	
	
	
	for( RS = IS = 0.0, m = 1; m <= M; m++)/* нулевое приближение */
	{
		RX[m] = RA[m][1+M];
		IX[m] = IA[m][1+M];

		RS += RX[m] * V[m];
		IS += IX[m] * V[m];
		
	}
		
	S = sqrt( RS * RS + IS * IS );

	for( i = 1; i <= It; i++)/* i-ое приближение */
	{
		S1 = S;
		for( m = 1; m <= M; m++)
		{
			Rx = RA[m][1+M];
			Ix = IA[m][1+M];
			for( n = 1; n <= M; n++)
			{
				Re = RA[m][n];
				Im = IA[m][n];

				BR = RX[n];
				BI = IX[n];

				Rx -= Re * BR - Im * BI;
				Ix -= Im * BR + Re * BI;
			}/* n */

			RX[m] = Rx;
			IX[m] = Ix;
		}/* m */

		for( RS = IS = 0.0, m = 1; m <= M; m++)
		{
			RS += RX[m] * V[m];
			IS += IX[m] * V[m];
		}

		S = sqrt( RS * RS + IS * IS );
	
	//	printf ("\n    n= %4d  %5.3e  %5.3e  %5.3e  %5.3e",i,RS,IS,S1,S);
	//	fprintf( out,"\n    n= %4d  %5.3e  %5.3e  %5.3e  %5.3e",i,RS,IS,S1,S);
		

		if ( fabs( S1 / S - 1.0 ) < ftol )
		{
			for( m = 1; m <= M; m++)
			{
				RA[m][1+M] = RX[m];
				IA[m][1+M] = IX[m];
			}
			break;
		}

		if ( i == It )
		{
			printf ("\n  Иттерации не сходятся");
		
			fmat_1d( (void*)RX, 1, sizeof( double) );
			fmat_1d( (void*)IX, 1, sizeof( double) );

		
		}

	}/* i */
	fmat_1d( (void*)RX, 1, sizeof( double) );
	fmat_1d( (void*)IX, 1, sizeof( double) );
}         /* end ITRSLAU */


