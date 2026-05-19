
#include "newlib.h"

void mult_a( double *x,double *y,double *z )//компл умножение
{
	z[1] = x[1] * y[1] - x[2] * y[2];
	z[2] = x[1] * y[2] + y[1] * x[2];
}

void div_a( double *x,double *y,double *z )//компл деление
{
	double zz;
	zz = 1.0 / ( y[1] * y[1] + y[2] * y[2] );
	z[1] = (x[1] * y[1] + x[2] * y[2]) * zz;
	z[2] = (-x[1] * y[2] + y[1] * x[2]) * zz;
}

void subt_a( double *x,double *y,double *z )//компл вычитание
{
	z[1] = x[1] - y[1];
	z[2] = x[2] - y[2];
}

void addit_a( double *x,double *y,double *z )//компл сложение
{
	z[1] = x[1] + y[1];
	z[2] = x[2] + y[2];
}

double mod( double *x,double *y)//модуль разности двух компл чисел
{
	double xx[3];
	xx[1] = x[1] - y[1];
	xx[2] = x[2] - y[2];
	return (sqrt(xx[1] * xx[1] + xx[2] * xx[2]));
}