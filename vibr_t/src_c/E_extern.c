
#include "VIBRr.h"

void E_extern ( double *kx_,double *ky_,double *kz_ )//падающее поле
{ // calc. external downward ElMagn-field
	// A.L. should extend for arbitrary illumination
	
	double k1,Z0_;
	//int P;

	k1 = sqrt(EPS_r[1]);
	tet0 = (fabs(TET) > 0.000001) ? TET_r : 1.0e-17;//”гол падени€

	teta = (H_[0] < 0.5) ? tet0 : asin( sin(tet0) / n_[1] );//”гол преломлени€

	cotet = cos(teta);
	sitet = sin(teta);

	cofi = cos(FI_r);
	sifi = sin(FI_r);
	
	kx = k1 * cofi * sitet;
	ky = k1 * sifi * sitet;
	kz = k1 * cotet;
	
	Z0_ = Z0 / n_[1];

	if (Pol)
	{ // P-polarization for Illum.azimuth along Y-axis

		Ex_0 = cofi * cotet;// ishodnoe pole
		Ey_0 = sifi * cotet;
		Ez_0 = -sitet;

		Hx_0 = sifi / Z0_;
		Hy_0 = -cofi / Z0_;
		Hz_0 =  0.0;
	}//Pol
	else
	{ // S-polarization for Illum.azimuth along Y-axis
		Ex_0 = sifi;// ishodnoe pole
		Ey_0 = -cofi;
		Ez_0 = 0.0;

		Hx_0 = - cofi * cotet / Z0_;
		Hy_0 = - sifi * cotet / Z0_;
		Hz_0 = sitet / Z0_;
	}//!Pol
	
	*kx_ = kx;
	*ky_ = ky;
	*kz_ = kz;
	
}//   E_extern

