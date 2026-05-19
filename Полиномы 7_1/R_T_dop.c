
#include "VIBRr.h"

//void Y(double gam2, int nu, double  *rYnu, double *iYnu);

double rYnu, iYnu ;

void R_T_dop(int PP, int m_u, int nu, double alf, double *rR_, double *iR_, double *rT_, double *iT_)
{
	double  gam, Q, Jmu, rR22, iR22,  C, rR1, iR1,  rR2, iR2,kap3,
		gam2, x, y, rR, iR, rTA, iTA, u, w;


	F_GRIN_dop(PP, alf, &rR, &iR, &rTA, &iTA);

	u = 1.0 - (rR * rR + iR * iR + rTA * rTA + iTA  * iTA);

	Jmu = J(m_u - 1, alf * L);

	gam2 = alf * alf - k2;
	gam = sqrt(-gam2);
	kap3 = sqrt(-alf * alf + k2 * eps[3] * mu__[3]);
	
	y = gam2 + alf * alf;
	x = (r_tau[1] + r_tau[2] * y) * Jmu;
	y = (i_tau[1] + i_tau[2] * y) * Jmu;
		
	
	Y(gam2, nu, &rYnu, &iYnu);
	mult(x,y, rYnu, iYnu,  &rR1, &iR1);
	//_________________

	Q = -alf * Q_mu(m_u - 1, L, -1, alf * L);
	
	rR2 = Q *  rYnu;
	iR2 = Q *  iYnu;//R2_1

	//_________________________
	C = 1.0 / (H * sqrt(Pi_2));
	

	Q = gam * Jmu * C;
			
	w = -2.0 *  gam * H + Pi05 * (double)(nu);
		
	y = sin(w);
	x = cos(w);

	w = - Pi05 * (double)(nu);
	y += sin(w);
	x += cos(w);


	//___________________________

	/*w = gam * H;

	u = cos(w);
	w = sin(w);

	mult(x, y, u, w, &x, &y);*/
	//_____________________________

	//mult(x, y, 0.0, 1.0, &x, &y);
	rR22 = x * Q;
	iR22 = y * Q;//R2_1

	mult(rR2 + rR22, iR2 + iR22, r_tau[2], i_tau[2], &x, &y);
		
	rR1 += x;
	iR1 += y;
	
	mult(rR, iR, rR1, iR1, &x, &y);
	*rR_ = x;
	*iR_ = y;
		
	mult(rTA, iTA, rR1, iR1, &x, &y);
	*rT_ = x;
	*iT_ = y;

		

}//R_T_dop

