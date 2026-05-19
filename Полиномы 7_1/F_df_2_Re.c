
#include "VIBRr.h"

void F_df_2_Re(int sp, //0-s-pol, 1-p-pol
	double *eps, double *mu, double *R, double *T)
{
	double  gam2,gam, kk, kap, kap2, kap_mu, C, S, rd, id, u, v;

	if (kbhit())
		if (getch() == ESCape)
		{
			fflush(out);
			exit(1);
		}
	
	
	kk = kx * kx;
	gam2 = k2 - kk;
	gam = sqrt(gam2);

	kap2 = k2 * eps[1]  * mu[1]- kk;
	kap = sqrt(kap2);
	
	kap_mu = (!sp) ? mu[1] : eps[1];
	kap_mu = kap / kap_mu;
	
	C = kap * h;
	S = sin(C);
	C = cos(C);

	rd = -(kap_mu * kap_mu + gam2) * S;
	id = 2.0 * kap_mu * gam * C;

	u = (kap_mu * kap_mu - gam2) * S;
	divide(u, 0.0, rd, id, &u, &v);
	*R = (u * u + v * v);

	v = 2.0 * kap_mu * gam;;
	divide(0.0,v, rd, id, &u, &v);
	*T = (u * u + v * v);
	
}//F_df
