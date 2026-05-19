// Функция нахождения нуля функции F_df
//  на интервале (COM, FN) с шагом ST.
// Возможно как  a<0,  так и a>0. Всегда  ST>0
//ik = 0 при входе. ik = 1 после нахождения перемены знака
//ftol - относительная точность

#include "newlib.h"

int sig(double x);
 
double GENERAL_L(double COM,double ST,double FN, double ftol,double (*F_df)(double),int *ik)
{
	int FIN;
	double fe,a,b,fa,fb,ff,s;
	
	s = (double)sig(FN - COM);
	FIN = 0;

	*ik = 0;
	a = COM;
	fa = F_df(a);

	met: b = a + ST * s;
	if ( b * s >= FN * s  )
	{
		b = FN;
		FIN = 1;
	}
	
	fb = F_df(b);
	
	if ( sig(fa) == sig(fb) )
	{
		a = b;
		fa = fb;
		if (!FIN)
			goto met;
		else
			return (0.0);
	}
	else
	{
		*ik = 1;
		fe = HORKUZ( a, b, fa, fb, F_df, &ff, ftol );
	}

	return (fe);


}   /*  GENERAL  */



int sig(double x)
{
	int y;
	y = ( x > 0.0 ) ? 1 : ( x < 0.0 ) ? -1 : 0;
	return (y);
}
