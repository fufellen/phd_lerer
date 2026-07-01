
#include "VIBRr.h"


void F_1( void ) 
{   
	double z,u1,v1,u,v;
	   	
	mult(EPS_r[1],EPS_i[1],EPS_r[2],EPS_i[2], &u,&v);
	u1 = EPS_r[1] + EPS_r[2];
	v1 = EPS_i[1] + EPS_i[2];

	divide(u,v, u1,v1, &u,&v);
	
	sq_compl_(u,v,&u,&v);//proverit
	

	z = (is[2]) ? sqrt(EPS_r[1]) : sqrt(EPS_r[2]);

	printf ("      %8.5lf   %8.5lf",u / z,-v/u);
	fprintf( out,"      %8.5lf   %8.5lf",u / z,-v/u);
	if (graf2)
		fprintf( out2,"    %8.5lf   %8.5lf",u ,-v/u);
	
}/* F */

