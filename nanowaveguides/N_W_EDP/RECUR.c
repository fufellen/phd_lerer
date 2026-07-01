
#include "VIBRr.h"


void RECURS(int N_d,int i,double **rte,double **ite,double **rqu,double **iqu, double *rG1,double *iG1, double *rG3,double *iG3)// 
{ 
	int N,n,bu,Nc;
	double u,v,rt[111],it[111],rq[111],iq[111],rq2[111],iq2[111],rA[111],iA[111],
		rD[111],iD[111],rfi,ifi,u1,v1,rP,iP;
	
	//здесь N_d -точка истока (в процедуре F_df точка истока - N_d1 ).

	for ( n=1; n <= Nsl + 1; n++ )  
	{	
		rt[n] = rte[i][n];
		it[n] = ite[i][n];
		rq[n] = rqu[i][n];
		iq[n] = iqu[i][n];

		mult(rq[n],iq[n], rq[n],iq[n], &rq2[n],&iq2[n]);
	}//n

	bu = (N_hol_m > 1) ? 1 : 0;
	
	// Nijnya FG- M3
		
	N = Nsl;
	rA[N+1] = iA[N+1] = iA[N] = 0.0;
	rA[N] = 1.0;
	for ( n=N; n > N_d; n-- )
	{
		rfi = rt[n] + rt[n+1];
		ifi = it[n] + it[n+1];
		
		mult(rA[n],iA[n], rfi,ifi, &u,&v);

		mult(rA[n+1],iA[n+1], rq2[n+1],iq2[n+1], &u1,&v1);

		rA[n-1] = u - u1;
		iA[n-1] = v - v1;
	}//n
		
	divide(1.0,0.0, rA[N_d],iA[N_d], &rP,&iP);//P = P * C

	for ( n=N_d + 1; n <= N; n++ )
	{	
		mult(rP,iP, rq[n],iq[n], &rP,&iP);
		mult(rA[n],iA[n], rP,iP, &rD[n],&iD[n]);
		
		if (n <= N_hol_m+1 && bu)
		{
			D_b[1][i][n][N_d] = rD[n];
			D_b[2][i][n][N_d] = iD[n];
		}
				
	}//n

	n = N_discon[N_hol] -1;
	//if (N_d <= N_hol && bu)
	if (N_d <= n && bu)
	{
		D_b[1][i][N_d][N_d] = 1.0;
		D_b[2][i][N_d][N_d] = 0.0;
	}
	

	mult(rD[N_d + 1],iD[N_d + 1], rq[N_d + 1],iq[N_d + 1], &u1,&v1);

	*rG3 = rt[N_d + 1] - u1;
	*iG3 = it[N_d + 1] - v1;

	DEN1[1] = rD[N];
	DEN1[2] = iD[N];
	// END Nijnya FG__________________________________________________________________________

	//// Verhnya FG- M1___________________________________________________________________________

	//N = N_di[N_d] - 1;
	
	N = N_d - 1;//здесь N_d -точка истока (в процедуре F_df точка истока - N_d1 ). 
	
	Nc = N_discon[1] - 1;
	
	//if (N == Nc)
	if (N == 1)
	{
		rD[0] = iD[0] = iD[1] = 0.0;
		rD[1] = 1.0;

		//rD[Nc-1] = iD[Nc-1] = iD[Nc] = 0.0;
		//rD[Nc] = 1.0;

		*rG1 = rte[i][N];
		*iG1 = ite[i][N];

	}//N == 1
	else/**/
	{			
		rA[0] = iA[0] = iA[1] = 0.0;
		rA[1] = 1.0;
		for ( n=1; n <= N - 1; n++ )
		{
			rfi = rt[n] + rt[n+1];//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			ifi = it[n] + it[n+1];
		
			mult(rA[n],iA[n], rfi,ifi, &u,&v);

			mult(rA[n-1],iA[n-1], rq2[n],iq2[n], &u1,&v1);

			rA[n+1] = u - u1;
			iA[n+1] = v - v1;
		}//n
				
		divide(1.0,0.0, rA[N],iA[N], &rP,&iP);//P = P * C

		for ( n=N-1; n >= 1; n-- )
		{	
			mult(rP,iP, rq[n+1],iq[n+1], &rP,&iP);

			mult(rA[n],iA[n], rP,iP, &rD[n],&iD[n]);
		
			D_t[1][i][n][N_d] = rD[n];
			D_t[2][i][n][N_d] = iD[n];
		}//n

		D_t[1][i][N][N_d] = 1.0;
		D_t[2][i][N][N_d] = 0.0;
		
		mult(rD[N-1],iD[N-1], rq[N],iq[N], &u1,&v1);
		*rG1 = rt[N] - u1;
		*iG1 = it[N] - v1;

	}//N != 1
		
	//DE1[1] = rD[Nc];//Ploxo dly D_n pri mz=20,21,22 Rehenie!
	//DE1[2] = iD[Nc];
	DE1[1] = rD[1];// znaki protivopoloj
	DE1[2] = iD[1];
	
	// END Verhnya FG
	
}//RECURS
		
				
			