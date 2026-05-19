
#include "VIBRr.h"

void  Imod_K(int n, double z, double *I, double *K)//Градштейн с. 980
//модифицированные ф. Бесселя с полуцелым индексом I_n+1/2(z) / z^1/2 и  K_n+1/2(z) / z^1/2
{
	double A,B,Bk, a, b, e,ep,s;
	int k,h=0;
	
	A = B = a = b = 1.0;
	s = -1.0;
	for (k = 0; k < n; k++)
	{
		s = -s;

		Bk = (double)((n + k + 1) *  (n - k)) / (double)(k + 1);
		Bk *= 0.5 / z;
		
		b *= Bk;
		a *= -Bk;

		A += a;
		B += b;
	}//k

	
	if (!h)//с экспонентами
	{
		e = exp(-z);
		ep = exp(z);

		*K = sqrt(0.5 * Pi) / z * e * B;

		*I = 1.0 / (sqrt(2.0 * Pi) * z) * (A * ep + e * s * B);
	}
	else// бкз экспонент
	{
		e = exp(-2.0 * z);
		
		*K = sqrt(0.5 * Pi) / z * B;

		*I = 1.0 / (sqrt(2.0 * Pi) * z) * (A + e * s * B);
	}
}//Imod_K

double  I_ex(int n, double z)//Градштейн с. 980
// модифицированная ф. Бесселя с полуцелым индексом J_n+1/2(z) / z^1/2 * exp(-z), n>=0
{
	double I,A, B, Bk, a, b, e, s;
	int k;

	A = B = a = b = 1.0;
	s = -1.0;
	for (k = 0; k < n; k++)
	{
		s = -s;

		Bk = (double)((n + k + 1) *  (n - k)) / (double)(k + 1);
		Bk *= 0.5 / z;

		b *= Bk;
		a *= -Bk;

		A += a;
		B += b;
	}//k

	e = exp(-2.0 * z);
	
	I = 1.0 / (sqrt(2.0 * Pi) * z) * (A + e * s * B);
	return I;
}//I_ex

double  J(int n, double z)//Градштейн с. 980
// ф. Бесселя с полуцелым индексоь J_n+1/2(z) / z^1/2, n>=0
{
	double A, Bk, a, b,  s;
	int k;

	b = 1.0;
	A = cos(z + 0.5 * Pi * (double)(-n - 1));
	
	for (k = 0; k < n; k++)
	{
		s = cos(z + 0.5 * Pi * (double)(-n +k));

		Bk = (double)((n + k + 1) *  (n - k)) / (double)(k + 1);
		Bk *= 0.5 / z;

		b *= Bk;
		a = b * s;

		A += a;
		
	}//k

	
	return (sqrt(2.0 / Pi) / z * A);
}//J

double  N(int n, double z)//Градштейн с. 980
// ф. Бесселя с полуцелым индексоь N_n+1/2(z) / z^1/2, n>=0
{
	double A,  Bk, a, b,  s;
	int k;

	A = b = 1.0;
	A = cos(z + 0.5 * Pi * (double)n);

	for (k = 0; k < n; k++)
	{
		s = cos(z + 0.5 * Pi * (double)(n + k + 1));

		Bk = (double)((n + k + 1) *  (n - k)) / (double)(k + 1);
		Bk *= 0.5 / z;

		b *= Bk;
		a = b * s;

		A += a;

	}//k


	return (sqrt(2.0 / Pi) / z * A * cos(Pi * (double)(n - 1)));
}//N

double  Pn(int n, double z)//Градштейн с. 1040 //Полиномы Лежандра
{
	double *P,k_;
	int k;

	P = (double*)amat_1d(0,n, sizeof(double));
	P[0] = 1.0;
	
	if (n)
	{
		P[1] = z;

		for (k = 1; k < n; k++)
		{
			k_ = (double)k;
			P[k + 1] = ( (2.0 * k_ + 1.0) * z * P[k] - k_ * P[k - 1]) / (k_ + 1.0);
		}//k
	}//n
	
	k_ = P[n];
	
	fmat_1d((void*)P, 0, sizeof(double));
	return k_;
}//Pn


double bess_dr_real(int n, double x, double eps, double nu, int i)//Функция Бесселя Jn+nu(x)/ x^nu

//i=0 - Модифицированные функция Бесселя
//i=1 - Функция Бесселя

{
	double term, sum, te, x2, n_;
	int s;

	if (fabs(x) <= 1.0e-7 && n >= 0)
	{
		if (!n)
			return(pow(0.5, nu) / gamma(1.0 + nu));
		else
			return(0.0);
	}
	else
	{
		n_ = (double)n;

		x2 = 0.25 * x * x;
		if (i)
			x2 = -x2;

		sum = term = pow(0.5, nu + n_) * pow(x, (double)n) / gamma(1.0 + nu + n_);
		s = 1;
		do
		{
			term *= x2 / (s * (s + n_ + nu));
			sum += term;
			te = fabs(term / sum);
			s++;
		} while (te > eps);
		return (sum);
	}
}  /* end bess */


double gamma(double x)
{
	double t, r, a[21];
	int i;

	a[0] = 1.0;
	a[1] = 0.577215664901532;
	a[2] = -0.655878071520253;
	a[3] = -0.042002635034095;
	a[4] = 0.166538611382291;
	a[5] = -0.042197734555544;
	a[6] = -0.009621971527876;
	a[7] = 0.007218943246663;
	a[8] = -0.001165167591859;
	a[9] = -0.000215241674114;
	a[10] = 0.000128050282388;
	a[11] = -0.000020134854780;
	a[12] = -0.000001250493482;
	a[13] = 0.000001133027231;
	a[14] = -0.205633841e-6;
	a[15] = 0.6116095e-8;
	a[16] = 0.5002007e-8;
	a[17] = -0.1181274e-8;
	a[18] = 0.104342e-9;
	a[19] = 0.7782e-11;
	a[20] = -0.3696e-11;

	t = 1.0 / x;
	if (x <= 1.0)
		goto m2;
	do
	{
		t *= x;
		x -= 1.0;
	} while (x > 1.0);

m2:
	r = a[20];
	for (i = 19; i >= 0; i--)
		r = a[i] + x * r;

	return(t / r);
}//gamma