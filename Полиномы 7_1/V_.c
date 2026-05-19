
#include "VIBRr.h"


void V_(int PP, int n, int nu, double alf, double *rV, double *iV)
{
	double  gam, gam2, x, y, u, w, gH;
		
	gam2 = alf * alf - k2;

	if (gam2 >= 0.0)
	{
		gam = sqrt(gam2);
		gH = gam * H;
		u = I_ex(n - 1, gH) * I_ex(nu - 1, gH) / gam;
		w = Pi05 * (double)(n - nu);
		x = cos(w) * u;
		y = sin(w) * u;
	}
	else
	{
		gam = sqrt(-gam2);
		gH = gam * H;
		u = J(n - 1, gH) * J(nu - 1, gH) / gam;

		w = -2.0 * gH + Pi * (double)(n - 1);

		x = sin(w) * u;
		y = -cos(w) * u;
	}
	
	*rV = x;
	*iV = y;

}