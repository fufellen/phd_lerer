
#include "VIBRr.h"


void compl_cos_sin(double x, double y, double *Rec, double *Imc, double *Res, double *Ims)
{
	double  s,sh,c,ch;
	
	s = sin(x);
	c = cos(x);
	
	sh = sinh(y);
	ch = cosh(y);

	*Rec = c * ch;
	*Imc = -s * sh;

	*Res = s * ch;
	*Ims = c * sh;
}//cos_sin

void compl_sqrt(double x, double y, double *Re, double *Im)
{
	double  mod,fi;

	mod = sqrt(x * x + y * y);
	mod = sqrt(mod);

	fi = 0.5 * atan2(y, x);

	*Re = mod * cos(fi);
	*Im = mod * sin(fi);
}//compl_sqrt

exp_compl(double x, double y, double *re, double *ie)//exp(iz)
{
	double r_ex;
	r_ex = exp(-y);
	*re = r_ex * cos(x);
	*ie = r_ex * sin(x);
}//exp_compl