
#include "VIBRr.h"


double F_(double *KZ)
{
	double F_C[3];

	if( kbhit() )
	   if( getch() == ESCape )
		 {
		 fflush( out );
		 exit( 1 );
		 }
		 
	cal++;

	F_Compl( KZ,F_C);
	
	//return (sqrt (F_C[1] * F_C[1] + F_C[2] * F_C[2]));
	//return (fabs (F_C[1]) * fabs(F_C[2]));
	return (fabs (F_C[1]* F_C[2]));
	
} 
