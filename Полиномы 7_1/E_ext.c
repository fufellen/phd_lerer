
#include "VIBRr.h"


void E_ext(int p, int mn_,double *rR, double *iR, double *rTA, double *iTA)
{
	double kx, gam, kap_1[4], kap[4],S,C,rA,iA,rP,iP,x,ret,iet,y,rd,id,mn;
	int m;
	
	kx = k * sin(TET);
	gam = k * cos(TET);

	for (m = 2; m <= 3; m++)
	{
		kap[m] = sqrt(k2 * epsG[m] * muG[m] - kx * kx);

		C = (p) ? epsG[m] : muG[m];
		kap_1[m] = kap[m] / C;
	}//m

	x = kap[2] * h2G;
	S = kap_1[2] / sin(x);
	C = cos(x);

	
	divide(S, 0.0, S * C, kap_1[3], &rA, &iA);

	mult(S, 0.0, -rA + C, -iA, &rP, &iP);

	//printf("\n s(anal)  S=%6.4lf   C=%6.4e    rA=%6.4e    iA=%6.4e     rP=%6.4e    iP=%6.4e ", S, C, rA, iA, rP, iP);
	
	ret = iP / gam;
	iet = -rP / gam;

	divide(1.0,0.0, 1.0 + ret,iet, &rd,&id);

	mult(1.0 - ret, -iet, rd,id, &x, &y);

	//printf("\n s(anal)  ret=%6.4lf   iet=%6.4e    rd=%6.4e    rd=%6.4e     rR=%6.4e    iR=%6.4e ", ret,iet,rd,id,x,y);
	*rR = x;
	*iR = y;
	
	mn = (mn_) ? 2.0 * sqrt(kap_1[3] / gam) : 2.0;
	mult(mn * rA, mn * iA, rd, id, &x, &y);

	*rTA = x;
	*iTA = y;

	
}//E_ext
