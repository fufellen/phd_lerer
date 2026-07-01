
#include "VIBRr.h"

void sq_compl(double x,double y,double *u,double *v)
//комплексный корень. 0 <= fi <= 2 * pi. !!!! y >= 0.  Для exp(-w*z), z>0 complex root. 0 <= fi <= 2 * pi. !! y> = 0. For exp (-w * z), z> 0
//исправил, справедливо для первой ветви корня, т.е. при 0 <= fi <= 2 * pi при любых y.
{
	double r,a;
	
	a = SQ_05;
	r = sqrt(x * x + y * y);
	*u = a * sqrt(r + x);
	
	//if (y < 0.0 )
		//*u = - *u;разрез вдоль положительной полуоси
		
	*v = a * sqrt(r - x);
	if (y < 0.0 )
		*v = - *v;//разрез вдоль отрицательной полуоси
}//sq_compl

void sq_compl_(double x,double y,double *u,double *v)
//комплексный корень. -pi <= fi <=  pi. y <= 0.
{
	double r,fi;
	
	r = sqrt(x * x + y * y);
	r = sqrt(r);
	fi = 0.5 * atan2(y,x);

	*u = r * cos(fi);
	*v =  r * sin(fi);
	/*a = SQ_05;
	r = sqrt(x * x + y * y);
	*u = a * sqrt(r + x);
	*v = -a * sqrt(r - x);
	*/
}//sq_compl

void SH_CTH(double x,double y,double *us,double *vs,double *uc,double *vc)
//us+i*vs = 1.0 / sh(x+iy)  uc+i*vc = cth(x+iy)
{
	double r,R,R2,c,s,c2,s2,zn;
	
	r = exp(-x);
	R = r * r;
	R2 = R * R;

	s = sin(y);
	c = cos(y);
	s2 = 2.0 * s * c;
	c2 = c * c - s * s;
	zn = 1.0 / (1.0 + R2 - 2.0 * R * c2);

	*uc = zn * (1.0 - R2);
	*vc = -zn * 2.0 * R * s2;

	zn *= 2.0 * r;
	*us = zn * (1.0 - R) * c;
	*vs = -zn * (1.0 + R) * s;
}

void C_sqrt(double x,double y, double *u,double *v)
{
	double r,fi;
	
	r =  sqrt(sqrt(x * x + y * y));
	fi = 0.5 * atan2(y,x);
	*u = r * cos(fi);
	*v = r * sin(fi);
	
}//C_sqrt

void C_exp(double x,double y, double *u,double *v)
{
	double r;
	
	r =  exp(x);
	*u = r * cos(y);
	*v = r * sin(y);	
}//C_exp

void C_sh(double x,double y, double *u,double *v)//0.5*(1-exp(-2*z)), z=x+i*y
{
	double r;
	if (x > 5.0)
	{
		*u = 0.5;
		*v = 0.0;
	}
	else
	{
		r =  0.5 * exp(- 2.0 * x);
		*u = 0.5 - r * cos(2.0 * y);
		*v = r * sin(2.0 * y);	
	}
}//C_sh

void C_ch(double x,double y, double *u,double *v)//0.5*(1+exp(-2*z)), z=x+i*y
{
	double r;
	if (x > 5.0)
	{
		*u = 0.5;
		*v = 0.0;
	}
	else
	{
		r =  0.5 * exp(- 2.0 * x);
		*u = 0.5  + r * cos(2.0 * y);
		*v = -r * sin(2.0 * y);
	}
}//C_ch

void C_sh_ch(double x,double y, double *us,double *vs, double *uc,double *vc)
{
	double r,rc,rs;
	
	r =  0.5 * exp(- 2.0 * x);
	rc = r * cos(2.0 * y);
	rs = r * sin(2.0 * y);

	*us = 0.5 - rc;
	*vs = rs;

	*uc = 0.5 + rc;
	*vc = -rs;

}//C_sh_ch

void Cn_sh_ch(double rg_,double ig_,double z,double h, double *us,double *vs, double *uc,double *vc)
{
	double x,y,u,v,u1,v1,e;
		
	x = -rg_ * (h - z);
	y = -ig_ * (h - z);
	e = 0.5 * exp(x);
	u = e * cos(y);
	v = e * sin(y);
	
	x = -rg_ * (h + z);
	y = -ig_ * (h + z);
	e = 0.5 * exp(x);
	u1 = e * cos(y);
	v1 = e * sin(y);

	*us = u - u1;
	*vs = v - v1;

	*uc = u + u1;
	*vc = v + v1;

}//C_sh

double Max(double x, double y)
{
	return( (x > y) ? x : y); 
}



void sin_cos_C(double x,double y,double *us,double *vs,double *uc,double *vc)
//us+i*vs = sin(x+iy)  uc+i*vc = cos(x+iy)
{
	double r,r1,c,s,sh,ch;
	
	r = 0.5 * exp(y);
	r1 = 0.25 / r;
	sh = r - r1; 
	ch = r + r1;
	
	s = sin(x);
	c = cos(x);

	*us = s * ch;
	*vs = c * sh;
	
	*uc = c * ch;
	*vc = -s * sh;
}//sin_cos_C


void sin_C(double x,double y,double *us,double *vs)
//us+i*vs = sin(x+iy)
{
	double r,r1;
	
	r = 0.5 * exp(y);
	r1 = 0.25 / r;
	
	*us = sin(x) * (r + r1);
	*vs = cos(x) * (r - r1);	
}//sin_C


void sinh_C(double x,double y,double *us,double *vs)
//us+i*vs = sh(x+iy)
{
	double r,r1;
	
	r = 0.5 * exp(x);
	r1 = 0.25 / r;
	
	*us = (r - r1) * cos(y);
	*vs = (r + r1) * sin(y);	
}//sinh_C

void cosh_C(double x,double y,double *us,double *vs)
//us+i*vs = sh(x+iy)
{
	double r,r1;
	
	r = 0.5 * exp(x);
	r1 = 0.25 / r;
	
	*us = (r + r1) * cos(y);
	*vs = (r - r1) * sin(y);	
}//sinh_C