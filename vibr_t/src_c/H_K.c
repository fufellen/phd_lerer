
#include "cylindr.h"

void H_K(double FI, double TETA)
{
	double coTET,siTET,coFI,siFI;

	siFI = sin( FI * Pi / 180.0 );
	coFI = cos( FI * Pi / 180.0 );

	siTET = sin( TETA * Pi / 180.0 );
	coTET = cos( TETA * Pi / 180.0 );

	kx = siTET * coFI;
	ky = siTET * siFI;
	kz = coTET;

	if (E)
	{
		H0x = - siFI;
		H0y = coFI;
		H0z = 0.0;		
	}
	else
	{
		H0x = coTET * coFI;
		H0y = coTET * siFI;
		H0z = - siTET;
	}

	
}   /*  H_K  */


	