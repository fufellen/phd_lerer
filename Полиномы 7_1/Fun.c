
#include "VIBRr.h"

double Q_mu(int mu, double l, int zn, double x)
{
	double a;
	
	a = (zn > 0) ? Pi05 * (double)mu : -Pi05 * (double)mu;
	a += x;
	//a = x + Pi05 * (double)mu;
	return (sin(a) / (l * sqrt(Pi05)));

}//Q_mu

void YSH(int nu, double gam, double *u, double*v)
{
	double gH,a;
	
	gH = gam * H;
	
	a = Q_mu(nu - 1, H, -1, gH) - gam * J(nu - 1, gH);
	*u = a * cos(gH);
	*v = -a * sin(gH);
}//YSH

void U_(int n, int nu, int cn, int cnu, double gam2, double *rA, double *iA)
{
	double mn, gam, I, K, z, ru,rU, iU,x,y,u;
	int n_, nu_,var = 1;
	if (cn != cnu)
		*rA = *iA = 0.0;
	else
	{
	
		if (n >= nu)
		{
			n_ = n;
			nu_ = nu;
		}
		else
		{
			n_ = nu;
			nu_ = n;
		}

		iU = 0.0;
		rU = (n == nu) ? (double)cn / (gam2 * H * Pi * (0.5 + (double)n)) : 0.0;
		
	

		if (gam2 >= 0.0)
		{
			gam = sqrt(gam2);
			z = gam * H;

			mn = cos(0.5 * Pi * (double)(n_ + nu_));
			//mn = 1.0;
			mn *= -2.0 / Pi;
			
			Imod_K(n_, z, &I, &K);
			ru = mn * I / gam;
		
			Imod_K(nu_, z, &I, &K);
			ru *= K;
			rU += ru;
			//iU = 0.0;
		}//gam2 >= 0.0
		else
		{
			gam = sqrt(-gam2);
			z = gam * H;
			//bess_dr_real( int n, double x, double eps, double nu, int i )
			if (var == 1)
			{
				x = J(n_, z);
				y = N(nu_, z);
			}
			else
			{
				x = bess_dr_real(nu_, z, 1.0e-6, 0.5, 1);
				y = N(nu_, z);
			}
			mn = -(double)cnu;//วัÝ !!!-

			u = (var == 1) ? J(n_, z) : J(n_, z)/*bess_dr_real(n_, z, 1.0e-6, 0.5, 1)*/;
			/*if ( fabs(J(n_, z) - bess_dr_real(n_, z, 1.0e-6, 0.5, 1)) > 1.0e-5 )
			{
				printf("\n  %2d  z=%6.4e    %6.4e   %6.4e", n_, z, J(n_, z), bess_dr_real(n_, z, 1.0e-6, 0.5, 1));
			}*/

			mn *= u / gam;//วัÝ !!!
			
			if (var == 1)
			{
				y = mn * J(nu_, z);
				x = mn * N(nu_, z);
			}
			else
			{
				y = mn * bess_dr_real(nu_, z, 1.0e-6, 0.5, 1);
				x = mn * N(nu_, z);
			}
			
			rU += x;
			iU = y;
			
			
		}//gam2 >< 0.0
		*rA = rU * (double)cn;
		*iA = iU * (double)cn;
	/*	*rA = rU;
		*iA = iU;*/
	}

}//U_
