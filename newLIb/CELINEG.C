
/*  пpоцедуpа pешения системы линейных уpавнений
 *    с комpлексными коэффициентами
 *      ( с одинаpной точностью )
 *
 *   Обpащение:  CELINEF( n, m, RA, IA );
 *               n - количество пpавых частей;
 *               m - поpядок системы;
 *               RA[1:m,1:m+n] - действительная часть коэффициентов;
 *               IA[1:m,1:m+n] - мнимая часть коэффициентов;
 */

#include "newlib.h"

void CELINEG(int M, int N, double** RA, double** IA)
{
	int i,j,k,MN;
	int l;
	double AI,AR,BI,BR,CI,CR,D;
	int ALB;
	MN = M+N;
	for( k = 1; k <= N; k++)
	{
		D = 0.0;
		for( i = k; i <= N; i++)
		{
			BR = fabs(RA[i][k]) + fabs(IA[i][k]);
			if (BR > D )
			{
				D = BR;l = i;
			}
		}
		if (D < 1.0e-18 )
		{
		  printf("\n celineg: STOP!  max( fabs( a[i,j] ) ) = 0.0 \n");
		  exit(1);
		 }
		AR = RA[l][k]; AI = IA[l][k];
		ALB =  ( fabs(AR) < fabs(AI) ) ? 1 : 0;
		if (ALB )
		{
			BR = AR;
			AR = AI; AI = BR;
		}    /*  ALB  */
		BR = AI/AR; AI = 1.0/( BR * AI + AR);
		AR = AI * BR;
		if ( ALB == 0 )
		{
			BR = AR; AR = AI; AI = BR;
		}
		RA[l][k] = RA[k][k]; IA[l][k] = IA[k][k];
		for( j=k+1; j <= MN; j++)
		{
			CR = RA[l][j]; CI = IA[l][j];
			if ( D < ( fabs(CR) + fabs(CI) ) * 1.0e-18 )
			{
				printf("\n CELINEG:STOP! \n");
			/*    fprintf(out,"\n CELINEG:STOP! \n");  */
				exit(1);
			}
			RA[l][j] = RA[k][j]; IA[l][j] = IA[k][j];
			BR =  CI * AI + CR * AR;
			BI = CI * AR - CR * AI;
			RA[k][j] = BR;
			IA[k][j] = BI;
			for( i = k+1; i <= N; i++)
			{
				CR = RA[i][k]; CI = IA[i][k];
				RA[i][j] = RA[i][j] - CR * BR + CI * BI;
				IA[i][j] = IA[i][j] - CR * BI - CI * BR;
			}   /* I */
		}   /* j */
	}   /*  k  */
	for( l = N+1; l <= MN; l++)
	for( i = N-1; i >= 1; i--)
	{
		AR = RA[i][l];AI = IA[i][l];
		for( j = i+1; j <= N; j++)
		{
			BR = RA[j][l];BI = IA[j][l];
			CR = RA[i][j];CI = IA[i][j];
			AR = AR-CR*BR+CI*BI;
			AI = AI-CR*BI-CI*BR;
		}   /*  j  */
		RA[i][l] = AR; IA[i][l] = AI;
	}
}         /* end CELINEF */


