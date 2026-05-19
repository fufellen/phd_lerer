
#include "VIBRr.h"


double rYnu, iYnu ;

void F_df_dop(int PP, int m, int m_u, int n, int nu, double alf, double *rA, double *iA)
{
	double  gam, Q, Jm, Jmu, rP, iP, rYn, iYn, zn,  C, rV, iV, co, si, rA1, iA1, rA2, iA2, Jm_mu,
		gam2,  x, y, rR, iR, rTA, iTA,w;
	
	
	F_GRIN_dop(PP, alf, &rR, &iR, &rTA, &iTA);

	Jm = J(m - 1, alf * L);
	Jmu = J(m_u - 1, alf * L);
	Jm_mu = Jm * Jmu;
	
	Q = Q_mu(m_u - 1, L, -1, alf * L);
	Q *= -alf * Jm;
	
	gam2 = alf * alf - k2;

	Y(gam2, nu, &rYnu, &iYnu);
	Y(gam2, n, &rYn, &iYn);
	zn = -cos(Pi * (double)n);
	rYn *= zn;
	iYn *= zn;
		
	mult(rYn, iYn,  rYnu, iYnu, &rV,&iV);
		
	co = Pi05 * (double)(nu - 1);
	si = sin(co);
	co = cos(co);

	C = 1.0 / (H * sqrt(Pi_2));
	if (gam2 >= 0.0)
	{
		gam = sqrt(gam2);
		
		rV /= gam;
		iV /= gam;
		
		w = exp(-2.0 * gam * H) ;
		rP = (co * w - co) * C;
		iP = (si * w + si) * C;
		
	}
	else
	{
		gam = sqrt(-gam2);
		
		divide(rV,iV, 0.0,gam,&rV, &iV);
				
		w = -2.0 *  gam * H + Pi * (double)(nu - 1);
		
		rP = (cos(w) - co) * C;
		iP = (sin(w) + si) * C;
	
	}//gam2 < 0.0)
	

	y = gam2 + alf * alf;
	x = r_tau[1] + r_tau[2] * y;
	y = i_tau[1] + i_tau[2] * y;

	mult(x, y, rV, iV, &x, &y);
	rA1 = Jm_mu * x;
	iA1 = Jm_mu * y;
	//______________________________

	rA2 = rV * Q;
	iA2 = iV * Q;

	mult(rYn, iYn, rP, iP, &x, &y);

	rA2 += x  * Jm_mu;
	iA2 += y * Jm_mu;

	mult(rA2, iA2, r_tau[2], i_tau[2], &rA2, &iA2);
	
	mult(rR, iR, rA1 + rA2,iA1 + iA2, &x, &y);

	
	*rA = x;
	*iA = y;

}//F_df_dop

void Y(double gam2, int nu, double *rYnu, double *iYnu)
{
	double gam,u,w;
	if (gam2 >= 0.0)
	{
		gam = sqrt(gam2);
		
		u = I_ex(nu - 1, gam * H);
		w = -Pi05 * (double)(nu - 1);
		*iYnu = sin(w) * u;
		*rYnu = cos(w) * u;
	}
	else
	{
		gam = sqrt(-gam2);
		
		u = J(nu - 1, gam * H);

		w = -gam * H;

		*iYnu = sin(w) * u;
		*rYnu = cos(w) * u;
	}
}//Y