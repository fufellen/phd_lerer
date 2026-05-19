
#include "VIBRr.h"

double gamma( double x);
void pow_compl( double x,double y,int n, double *u,double *v );

//Вычисления функций Бесселя действительного и комплексного аргумента разложением в ряд

 double bess_dr_real( int n, double x, double eps, double nu, int i )//Функция Бесселя Jn+nu(x)/ x^nu
 
 //i=0 - Модифицированные функция Бесселя
 //i=1 - Функция Бесселя

 {
	double term, sum, te,x2,n_;
	int s;

	if ( fabs(x) <= 1.0e-7 && n >= 0)
	{
		if ( !n )
			return( pow( 0.5, nu ) / gamma( 1.0 + nu ) );
         else
			return( 0.0 );
	}
	else
	{
		n_ = (double)n;

		x2 = 0.25 * x * x;
		if (i)
			x2 = - x2;
		
		sum = term = pow( 0.5, nu+n_ ) * pow( x,( double )n ) / gamma( 1.0 + nu+n_ );
		s = 1;
		do
		{
			term *= x2 / ( s * ( s + n_ + nu ) );
			sum += term;
			te = fabs( term / sum );
			s++;
		} while ( te > eps );
		return ( sum );
	}
 }  /* end bess */

 void bess_dr_compl( int n, double x, double y, double eps, double nu, double *Re_j, double *Im_j )//Функция Бесселя Jn+nu(z)/ z^nu
 //z=x+iy
 //i=0 - Модифицированные функция Бесселя
 //i=1 - Функция Бесселя

 {
	double term[3], sum[3], te,x2,n_, m_z2, m_z,u,v,z2[3],k;
	int s;

	m_z2 = x * x + y * y;
	m_z = sqrt(m_z2);
	if ( m_z <= 1.0e-7 && n >= 0)
	{
		*Im_j = 0.0;
		
		*Re_j = ( !n ) ? pow( 0.5, nu ) / gamma( 1.0 + nu ) : 0.0;
		
	}
	else
	{
		n_ = (double)n;

		z2[1] = -0.25 * (x * x - y * y);
		z2[2] = -0.5 * x * y;
		
		
		te = pow( 0.5, nu+n_ )  / gamma( 1.0 + nu+n_ );

		pow_compl( x,y,n, &u,&v );
		sum[1] = term[1] = te * u;
		sum[2] = term[2] = te * v;
		
		s = 1;
		do
		{
			k =1.0 / ( s * ( s + n_ + nu ) );
			mult_a(term, z2, term);
			term[1] *= k;
			term[2] *= k;
			
			sum[1] += term[1];
			sum[2] += term[2];
			te = (fabs( sum[1]) > 1.0e-8) ? fabs( term[1] / sum[1] ) : 0.0;
			k = (fabs( sum[2]) > 1.0e-8) ? fabs( term[2] / sum[2] ) : 0.0;
			s++;
		} while ( te > eps ||  k > eps);
		
		*Re_j = sum[1];
		*Im_j = sum[2];
	}
 }  /* end bess */

double gamma( double x)
{
	double t,r,a[21];
	int i;

	a[0] = 1.0;
	a[1] =  0.577215664901532;
	a[2] = -0.655878071520253;
	a[3] = -0.042002635034095;
	a[4] =  0.166538611382291;
	a[5] = -0.042197734555544;
	a[6] = -0.009621971527876;
	a[7] =  0.007218943246663;
	a[8] = -0.001165167591859;
	a[9] = -0.000215241674114;
	a[10] =  0.000128050282388;
	a[11] = -0.000020134854780;
	a[12] = -0.000001250493482;
	a[13] =  0.000001133027231;
	a[14] = -0.205633841e-6;
	a[15] =  0.6116095e-8;
	a[16] =  0.5002007e-8;
	a[17] = -0.1181274e-8;
	a[18] =  0.104342e-9;
	a[19] =  0.7782e-11;
	a[20] = -0.3696e-11;

	t = 1.0 / x;
	if ( x <= 1.0 )
		goto m2;
	do
	{
		t *= x;
		x -= 1.0;
	} while ( x > 1.0 );

	m2:
	r = a[20];
	for ( i = 19; i >= 0; i-- )
		r = a[i] + x * r;

	return( t / r );
}//gamma

void pow_compl( double x,double y,int n, double *u,double *v )//z^n
{
	int i;
	double U,V;
	
	U = 1.0;
	V = 0.0;

	for ( i = 1; i <= abs(n); i++ )
		mult(U,V, x,y, &U,&V);
	
	if (n < 0)
		divide(1.0,0.0, U,V, &U, &V);

	*u = U;
	*v = V;
	
}//pow_compl