
void mult( double x1,double y1,double x2,double y2,double *u,double *v )
{
	*u = x1 * x2 - y1 * y2;
	*v = x1 * y2 + y1 * x2;
}

void divide(double x,double y,double x1,double y1,double *u,double *v)
{
	double z;
	z = 1.0 / ( x1 * x1 + y1 * y1 );
	*u = ( x * x1 + y * y1 ) * z;
	*v = ( -x * y1 + y * x1 ) * z;
}
