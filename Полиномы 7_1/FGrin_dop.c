
#include "VIBRr.h"


void F_GRIN_dop(int p, double alf, double *rR, double *iR, double *rTA, double *iTA)
{
	double gam, kap_1[4], kap[4], S, C, rA, iA, rP, iP, x, ret, iet, y, rd, id,  fin = 120.0;
	int m, bu[4];

	eps[1] = mu__[1] = 1.0;
	for (m = 1; m <= 3; m++)
	{
		kap[m] = -k2 * eps[m] * mu__[m] + alf * alf;
		if (kap[m] > 0.0)
		{
			bu[m] = 1;
			kap[m] = sqrt(kap[m]);
		}
		else
		{
			bu[m] = 0;
			kap[m] = sqrt(-kap[m]);
		}

		C = (p) ? eps[m] : mu__[m];
		kap_1[m] = kap[m] / C;
	}//m

	gam = kap[1];

	x = kap[2] * h2;
	if (bu[2] && x > fin)
	{
		ret = kap_1[2] / kap_1[1];
		*rR = (1.0 - ret) / (1.0 + ret);
		*iR = *rTA = *iTA = 0.0;
		return;
	}//bu[m] && x > fin

	if (bu[2])
	{
		S = kap_1[2] / sinh(x);
		C = cosh(x);
	}
	else
	{
		S = kap_1[2] / sin(x);
		C = cos(x);
	}

	if (bu[3])
	{
		rA = S / (kap_1[3] + S * C);
		rP = S * (-rA + C);
		iA = iP = 0.0;
	}
	else
	{
		divide(S, 0.0, S * C, kap_1[3], &rA, &iA);
		mult(S, 0.0, -rA + C, -iA, &rP, &iP);
	}
	
	if (bu[1])
	{
		ret = rP / kap_1[1];
		iet = iP / kap_1[1];
	}
	else
	{
		ret = iP / gam;
		iet = -rP / gam;
	}
		
	divide(1.0,0.0, 1.0 + ret,iet, &rd,&id);

	mult(1.0 - ret, -iet, rd, id, &x, &y);

	
	*rR = x;
	*iR = y;
	
	mult(2.0 * rA, 2.0 * iA, rd, id, &x, &y);

	
	*rTA = x;
	*iTA = y;

	
}//E_ext
