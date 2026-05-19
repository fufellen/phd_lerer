
#include "newlib.h"
//Поиск комплексных нулей методом Ньютона
//st << 1 - параметр для численного нахождения производной
//x[1] + i * x[2] - начальное приближение, i -мнимая единица. Результат записывается в x[]
//ff[1] + i * ff[2]- значение функции в нуле
//tol - относительная тчность
void NEW_C( double *x, double st, void func(double *z,double *f), double* ff, double tol)
{
	double b[3],x1[3],fma[3],fx[3],fb[3],bx[3],fbx[3],z[3],fpr[3];
	int n = 0;

	func(x,fx);

	fma[1] = fx[1];
	fma[2] = fx[2];
	
	do
	{
		//b[1] = x[1] + st;
		//b[2] = x[2] + st;
		b[1] = x[1] * (1.0 + st);
		b[2] = x[2] * (1.0 + st);
	
		func(b,fb);
		

		subt_a(b,x,bx);
		subt_a(fb,fx,fbx);
		div_a(fbx,bx,fpr);

		div_a(fx,fpr,z);

		subt_a(x,z,x1);
       		
		x[1] = x1[1];
		x[2] = x1[2];
		func(x,fx);
		
		n++;

	}while (sqrt(z[1] * z[1] + z[2] * z[2]) > tol && n < 100);

	
	ff[1] = fx[1] / fma[1];
	ff[2] = fx[2] / fma[2];
}

