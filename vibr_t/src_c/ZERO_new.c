
#include "VIBRr.h"

int sig(double x);
double GENERAL(double  COM,double step,double ftol,double (*F_df)(double) );
double F_m (double x);
double F_m_C (double x);

double Rnu,Rnu1;
int m;

void ZERO_B(void) 
{
	int by,n,i,nl,l,m_r;
	double A,x,COM,st,step,ftol=1.0e-8,J_;

	for ( n=1;n <= N_hol;n++ )
	{
		for ( nl = 1; nl <= N_l[n]; nl++ )
		{			
			Rnu = r_l[n][nl]; 
			Rnu1 = (nl == 1) ? 0.0 :r_l[n][nl-1];
			st = Rnu - Rnu1;
			step = 0.01 * Pi / st;

			m_r = (nl == 1) ? M_r_1 : M_r_;
		
			for ( m=0;m <= M_fi - 1;m++ )
			{
				COM = 0.1 / st;
				for (  l=1; l <= m_r; l++ )
				{
					i = (nl == 1) ? l : l + (nl - 1) * M_r_  + 1;
					
					N_lin[n][i] = nl;
					
					if (m + l == 1)
					{
						by = 1;
						x = 0.0;
					}
					else
					{
						by = 0;
						x = (nl == 1) ? GENERAL(COM,step,ftol,F_m) : GENERAL(COM,step,ftol,F_m_C);
					}
					
					j[m][i][n] = x;

					J[m][i][n] = J_ = jn(m,x * Rnu);
					
					if (i > M_r_1)
					{
						if (by)
						{
							A_k[m][i][n] = 0.0;
							J[m][i][n] = Y1[m][i][n] = 1.0;
						}
						else
						{
							A_k[m][i][n] = A = - ( jn(m-1,x * Rnu) - jn(m+1,x * Rnu) ) / ( yn(m-1,x * Rnu) - yn(m+1,x * Rnu) );
							J[m][i][n] += A * yn(m,x * Rnu);
							Y1[m][i][n] = jn(m,x * Rnu1) + A * yn(m,x * Rnu1);
						}
					
					}//i > M_r1

					COM = x + step;
			
					
				} /* nl */
			}//m
		}//nl
	}//n
	

} //ZERO_B



double GENERAL(double  COM,double step,double ftol,double (*F_df)(double) )
{
	double fe,a,b,fa,fb,ff;

	a = COM/* - Abet[1]*/;
	
	fa = F_df(a);

	met: b = a + step;
		
	fb = F_df(b);
	
	if ( sig(fa) == sig(fb) )
	{
		a = b;
		fa = fb;
		goto met;
	}
	else
	fe = HORKUZ( a, b, fa, fb, F_df, &ff, ftol );//v library Newlib
	
	return (fe);


}   /*  GENERAL  */



int sig(double x)
{
	int y;
	y = ( x > 0.0 ) ? 1 : ( x < 0.0 ) ? -1 : 0;
	return (y);
}



double F_m (double x)
{
	double  Fm;

	x *= Rnu;
	Fm = jn(m-1,x) -  jn(m+1,x);

	return (Fm);
} /* F_m */


double F_m_C (double x)
{
	double  Fm,y,Jp,Jp1,Np,Np1;

	y = x * Rnu1;
	x *= Rnu;

	Jp = jn(m-1,x) - jn(m+1,x);
	Np = yn(m-1,x) - yn(m+1,x);

	Jp1 = jn(m-1,y) - jn(m+1,y);
	Np1 = yn(m-1,y) - yn(m+1,y);

	Fm = Jp * Np1 - Jp1 * Np;

	return (Fm);
} /* F_m */

/* Proverca
y = jn(n-1,x) -  jn(n+1,x);
						
printf ("\n   %3d    %3d    	%12.9e	  %12.9e", n,i,x, y);
fprintf( out,"\n   %3d    %3d    	%12.9e	  %12.9e", n,i,x, y);
			
if ( graf1  )
	fprintf( out1,"\n   j[%1d][%2d] = %12.9e;", n,i,x);
Proverca end
*/
			

//printf ("\n   %3d    %3d        %3d    %3d    %3d", n,nl,m,l,i);
//printf ("\n   %3d    %3d    %3d    	%12.9e	  %12.9e	  %12.9e",m,i,n,J[m][i][n],Y[m][i][n], Y1[m][i][n]);
					//printf ("\n  %2d  %2d   %2d   %2d    	%5.3e	  %5.3e	  %5.3e",m,nl,i,n,J[m][i][n],Y[m][i][n], Y1[m][i][n]);
					