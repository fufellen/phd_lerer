
#include "VIBRr.h"


void RECURS_ZERO(int N_d,int i,double **rte,double **ite,double **rqu,double **iqu, double *rG1,double *iG1, double *rG3,double *iG3)// 
{ 
	int N,n,bu,Nsl_,N_d_;
	double u,v,rt[177],it[177],rq[177],iq[177],rq2[177],iq2[177],rA[177],iA[177],
		rD[177],iD[177],rfi,ifi,u1,v1,rP,iP,uu;
	//double db[3][177],dt[3][177];
	
	//здесь N_d -точка истока (в процедуре F_df точка истока - N_d1 ).

	for ( n=1; n <= Nsl + 1; n++ )  
	{	
		rt[n] = rte[i][n];
		it[n] = ite[i][n];
		rq[n] = rqu[i][n];
		iq[n] = iqu[i][n];

		mult_(rq[n],iq[n], rq[n],iq[n], &rq2[n],&iq2[n]);
	}//n

	bu = (N_hol > 1) ? 1 : 0;
	
	// Nijnya FG- M3

	Nsl_ = Nsl;
	/**/
	for ( n=N_d+1; n <= Nsl; n++ )
	{
		if (fabs(rq[n]) < 1.0e-20)
		{
			Nsl_ = n - 1;
			break;
		}
	}//n
	
	for ( n=Nsl+1; n >= 1; n-- )
		rA[n] = iA[n] = 0.0;
	
		
	N = Nsl_;
	rA[N+1] = iA[N+1] = iA[N] = 0.0;// zdes A=D v 1-oi recurentnoi formule
	rA[N] = 1.0;

	for ( n=N; n > N_d; n-- )
	{
		rfi = rt[n] + rt[n+1];
		ifi = it[n] + it[n+1];
		
		mult_(rA[n],iA[n], rfi,ifi, &u,&v);

		mult_(rA[n+1],iA[n+1], rq[n+1],iq[n+1], &u1,&v1);
		divide_(u - u1,v - v1, rq[n],iq[n], &u1,&v1);

		rA[n-1] = u1;
		iA[n-1] = v1;
	}//n
		
	uu = 1.0 / ( rA[N_d] *  rA[N_d] + iA[N_d] * iA[N_d]);
	rP = uu * rA[N_d];
	iP = -uu * iA[N_d];

	for ( n=N_d + 1; n <= N; n++ )
	{	
		mult_(rA[n],iA[n], rP,iP, &rD[n],&iD[n]);
		
		if (n <= N_hol+1 && bu)
		{
			D_b[1][i][n][N_d] = rD[n];
			D_b[2][i][n][N_d] = iD[n];
		}
		
	}//n

	if (N_d <= N_hol && bu)
	{
		D_b[1][i][N_d][N_d] = 1.0;
		D_b[2][i][N_d][N_d] = 0.0;
	}
	
	
	mult_(rD[N_d + 1],iD[N_d + 1], rq[N_d + 1],iq[N_d + 1], &u1,&v1);

	*rG3 = rt[N_d + 1] - u1;
	*iG3 = it[N_d + 1] - v1;
	DEN1[1] = rD[N];
	DEN1[2] = iD[N];
	// END Nijnya FG__________________________________________________________________________

	//// Verhnya FG- M1___________________________________________________________________________

	N = N_d - 1;//здесь N_d -точка истока (в процедуре F_df точка истока - N_d1 ). 

	if (N == 1)
	{
		rD[0] = iD[0] = iD[1] = 0.0;
		rD[1] = 1.0;

		*rG1 = rte[i][N];
		*iG1 = ite[i][N];

	}//N == 1
	else/**/
	{	
		for ( n=1; n <= N_d - 1; n++ )
		{
			if (fabs(rq[n+1]) < 1.0e-20)
			{
				N_d_ = n-1;
				break;
			}
			
		}//n

		N = N_d - 1;

		rA[0] = iA[0] = iA[1] = 0.0;
		rA[1] = 1.0;
		for ( n=1; n <= N - 1; n++ )
		{
			rfi = rt[n] + rt[n+1];
			ifi = it[n] + it[n+1];
		
			mult_(rA[n],iA[n], rfi,ifi, &u,&v);

			mult_(rA[n-1],iA[n-1], rq[n],iq[n], &u1,&v1);
			divide_(u - u1,v - v1, rq[n+1],iq[n+1], &u1,&v1);

			rA[n+1] = u1;
			iA[n+1] = v1;
			
		}//n
				
		

		uu = 1.0 / ( rA[N] *  rA[N] + iA[N] * iA[N]);
		rP = uu * rA[N];
		iP = -uu * iA[N];

		for ( n=N-1; n >= 1; n-- )
		{	
			mult_(rA[n],iA[n], rP,iP, &rD[n],&iD[n]);
		
			D_t[1][i][n][N_d_] = rD[n];
			D_t[2][i][n][N_d_] = iD[n];
		}//n

	
		D_t[1][i][N][N_d_] = 1.0;
		D_t[2][i][N][N_d_] = 0.0;
		
		mult_(rD[N-1],iD[N-1], rq[N],iq[N], &u1,&v1);
		*rG1 = rt[N] - u1;
		*iG1 = it[N] - v1;

	}//N != 1
		
	DE1[1] = rD[1];
	DE1[2] = iD[1];
	
	// END Verhnya FG
	
	
}//RECURS
		
				
			