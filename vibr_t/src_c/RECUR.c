
#include "VIBRr.h"


void RECURS(double Ro,int N_d,int i,double **rte,double **ite,double **rqu,double **iqu, double *rG1,double *iG1, double *rG3,double *iG3)// 
{ 
	int N,n,bu;
	double u, v, rt[177], it[177], rq[177], iq[177], rA[177], iA[177], 
		rpci[111], ipci[111], Hn[111],  ex, Ro_,
		rD[177], iD[177], rfi, ifi, u1, v1, u2, v2, rP, iP, uu;
	//double db[3][177],dt[3][177];
	
	//здесь N_d -точка истока (в процедуре F_df точка истока - N_d1 ).
	Ro_ = sqrt(Ro);

	for ( n=1; n <= Nsl + 1; n++ )  
	{	
		rt[n] = rte[i][n];
		it[n] = ite[i][n];
		rq[n] = rqu[i][n];
		iq[n] = iqu[i][n];
	}//n
	
	bu = (N_hol > 1) ? 1 : 0;
	
	// Nijnya FG- M3
		
	N = Nsl;
	
	rpci[N] = 1.0;
	rpci[N+1] = ipci[N+1] = ipci[N] = 0.0;
	Hn[N] = Hn[N+1] = 0.0;

	for (n = N; n >= N_d; n--)
	{
		Hn[n - 1] = Hn[n] + H_[n-1];

		divide_(dze[i][n+1], dze_im[i][n+1], dze[i][n], dze_im[i][n], &u1, &v1);
		mult_(rpci[n], ipci[n], u1 + 1.0, v1, &u, &v);
		rpci[n - 1] = 0.5 *u;
		ipci[n - 1] = 0.5 * v;
	}//n
	
	rA[N+1] = iA[N+1] = iA[N] = 0.0;// zdes A=D v 1-oi recurentnoi formule
	rA[N] = 1.0;

	
	for ( n=N; n > N_d; n-- )
	{
		rfi = rt[n] + rt[n+1];
		ifi = it[n] + it[n+1];
		
		mult_(rA[n],iA[n], rfi,ifi, &u,&v);//A[n]*fi
		/**/
		divide_(rpci[n], ipci[n],rpci[n - 1], ipci[n - 1], &u1, &v1);
		
		mult_(u,v, u1,v1, &u,&v);
		
		mult_(rA[n+1],iA[n+1], rq[n+1],iq[n+1], &u1,&v1);//A[n+1]*q

		/**/
		divide_(rpci[n+1], ipci[n+1], rpci[n - 1], ipci[n - 1], &u2, &v2);
		ex = exp(-2.0 * Ro_ * H_[n]);
		u2 *= ex;
		v2 *= ex;
		mult_(u1, v1, u2, v2, &u1, &v1);
		

		divide_(u - u1,v - v1, rq[n],iq[n], &u1,&v1);

		rA[n-1] = u1;
		iA[n-1] = v1;
	}//n
		
	mult_(rA[N_d], iA[N_d], rpci[N_d], ipci[N_d], &u1, &v1);
	uu = 1.0 / (u1 * u1 + v1 * v1);

	rP = uu * u1;
	iP = -uu * v1;

	for ( n=N_d + 1; n <= N; n++ )
	{	
		mult_(rP, iP, rpci[n], ipci[n], &u,&v);
		ex = exp(Ro_ * (Hn[n] - Hn[N_d]));
		u *= ex;
		v *= ex;
		
		mult_(rA[n],iA[n], u,v, &rD[n],&iD[n]);
		
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
	
	ex = exp(-Ro_ * H_[N_d]);
	mult_(rD[N_d + 1],iD[N_d + 1], ex * rq[N_d + 1],ex * iq[N_d + 1], &u1,&v1);

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
		rpci[1] = 1.0;
		rpci[0] = ipci[1] = ipci[0] = 0.0;
		Hn[0] = Hn[1] = 0.0;

		for (n = 1; n <= N; n++)
		{
			Hn[n + 1] = Hn[n] + H_[n];

			divide_(dze[i][n], dze_im[i][n], dze[i][n + 1], dze_im[i][n + 1], &u1, &v1);
			mult_(rpci[n], ipci[n], u1 + 1.0, v1, &u, &v);
			rpci[n + 1] = 0.5 * u;
			ipci[n + 1] = 0.5 * v;
		}//n


		rA[0] = iA[0] = iA[1] = 0.0;
		rA[1] = 1.0;

		for ( n=1; n <= N - 1; n++ )
		{
			rfi = rt[n] + rt[n+1];
			ifi = it[n] + it[n+1];
		
			mult_(rA[n],iA[n], rfi,ifi, &u,&v);

			divide_(rpci[n], ipci[n], rpci[n + 1], ipci[n + 1], &u1, &v1);
						
			mult_(u, v, u1, v1, &u, &v);

			mult_(rA[n-1],iA[n-1], rq[n],iq[n], &u1,&v1);

			divide_(rpci[n - 1], ipci[n - 1], rpci[n + 1], ipci[n + 1], &u2, &v2);
			ex = exp(-2.0 * Ro_ * H_[n - 1]);
			u2 *= ex;
			v2 *= ex;
			mult_(u1, v1, u2, v2, &u1, &v1);

			divide_(u - u1,v - v1, rq[n+1],iq[n+1], &u1,&v1);

			rA[n+1] = u1;
			iA[n+1] = v1;
			
		}//n
				
				
		mult_(rA[N], iA[N], rpci[N], ipci[N], &u1, &v1);
		uu = 1.0 / (u1 * u1 + v1 * v1);

		rP = uu * u1;
		iP = -uu * v1;


		for ( n=N-1; n >= 1; n-- )
		{	
			mult_(rP, iP, rpci[n], ipci[n], &u, &v);
			ex = exp(Ro_ * (Hn[n] - Hn[N]));
			u *= ex;
			v *= ex;
		
			mult_(rA[n],iA[n], u,v, &rD[n],&iD[n]);
			D_t[1][i][n][N_d] = rD[n];
			D_t[2][i][n][N_d] = iD[n];
		}//n

		D_t[1][i][N][N_d] = 1.0;
		D_t[2][i][N][N_d] = 0.0;
		
		ex = exp(-Ro_ * H_[N-1]);
		mult_(rD[N-1],iD[N-1], ex * rq[N],ex * iq[N], &u1,&v1);
		*rG1 = rt[N] - u1;
		*iG1 = it[N] - v1;

	}//N != 1
		
	DE1[1] = rD[1];
	DE1[2] = iD[1];
	
	// END Verhnya FG
	
	
}//RECURS
		
				
			