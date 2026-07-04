
#include "VIBRr.h"


 void CONST(int N_d,int i,double Ro, double *te,double *qu,double *fi1,double *fi3,double *del,double *S1,double *S3)
 {
	int n;
	double dze[3],gam[3],gam_[3],u,v,u1,v1,u2,v2,x,y,**rq,**iq,**rte,**ite,Dz[3],rga1,rga3,iga1,iga3;

	rte = ( double**)amat_2d( 1,3, 1,Nsl+1, sizeof( double) );
	ite = ( double**)amat_2d( 1,3, 1,Nsl+1, sizeof( double) );
	rq = ( double**)amat_2d( 1,3, 1,Nsl+1, sizeof( double) );
	iq = ( double**)amat_2d( 1,3, 1,Nsl+1, sizeof( double) );

	for ( n=1; n <= Nsl + 1; n++ )//Вачисление gamma[1:3] с черточкой и q, t Calculation of gamma [1:3] with a dash and qб t
	{
		u = Ro - k2_r[n];
		v = -k2_i[n];
		sq_compl(u,v, &gam[1],&gam[2]);
				
		if ( i == 1)
		{
			dze[1] = 1.0;
			dze[2] = 0.0;
		}
		else
			{
				if (i == 3)
					divide_(1.0,0.0, EPS_r[n],EPS_i[n], &x,&y);
				else
					divide_(EPS_r[n],EPS_i[n], u,v, &x,&y);
			
				dze[1] = x;
				dze[2] = y;
			}

		mult_(gam[1],gam[2], dze[1],dze[2], &gam_[1],&gam_[2]);
		
		if (n == 1 || n == Nsl + 1)
		{
			rte[i][n] = gam_[1];
			ite[i][n] = gam_[2];
			rq[i][n] = iq[i][n] = 0.0;
		}
		else
		{			
			u = gam[1] * H_[n-1];
			v = gam[2] * H_[n-1];
			SH_CTH(u,v,&u1,&v1,&u2,&v2);

			mult_(gam_[1],gam_[2], u1,v1, &u,&v);
			rq[i][n] = u;
			iq[i][n] = v;

			mult_(gam_[1],gam_[2], u2,v2, &u,&v);
			rte[i][n] = u;
			ite[i][n] = v;

		}

		if (n == N_d)
		{
			te[1] = rte[i][n];
			te[2] = ite[i][n];

			qu[1] = rq[i][n];
			qu[2] = iq[i][n];

			Dz[1] = dze[1];
			Dz[2] = dze[2];
		}
	}//n

	//RECURS(N_d,i,rte,ite,rq,iq, &rga1,&iga1, &rga3,&iga3);
	RECURS_ZERO(Ro,N_d,i,rte,ite,rq,iq, &rga1,&iga1, &rga3,&iga3);
	//________________________________________________
	fi1[1] = rga1 + te[1];
	fi1[2] = iga1 + te[2];

	fi3[1] = rga3 + te[1];
	fi3[2] = iga3 + te[2];
	
	//_________________________________


	mult_(fi1[1],fi1[2], fi3[1],fi3[2], &x,&y);
				
	u = rq[i][N_d] * rq[i][N_d] - iq[i][N_d] * iq[i][N_d];
	v = 2.0 * rq[i][N_d] * iq[i][N_d];//q2^2

	//del[1] = x - u;//old
	//del[2] = y - v;
	
	divide_(Dz[1],Dz[2], x - u,y - v, &del[1],&del[2]);//news
	//_________________________________
	
	S3[1] = DEN1[1];
	S3[2] = DEN1[2];
	S1[1] = DE1[1];
	S1[2] = DE1[2];

	fmat_2d( (void**)rte,1,3, 1, sizeof( double) );
	fmat_2d( (void**)ite,1,3, 1, sizeof( double) );
	fmat_2d( (void**)rq,1,3, 1, sizeof( double) );
	fmat_2d( (void**)iq,1,3, 1, sizeof( double) );


}//CONST