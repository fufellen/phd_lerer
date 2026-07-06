
#include "VIBRr.h"

/*
Ѕиблиотека материалов
¬ход Ц нм, кроме рентгена в eV.
¬ыход Ц действительна€ и мнима€ часть диэлектрической проницаемости.

is[n] = 1 медь
is[n] = 2 - серебро
is[n] = 3 - золото
is[n] = 4 - ZnO
is[n] = 5 - Si
is[n] = 6 - SiO2
is[n] = 7 - HfO
is[n] = 8 - AZO - база 1
is[n] = 9 - рентген база в eV  реализовано не во всех программах
is[n] = 10 AZO_2 - база 2
is[n] = 11 - AZO_wt_0_8_
is[n] = 12 - AZO_wt_2_0_
is[n] = 13 - AZO_wt_3_2_
is[n] = 14 - GZO_
is[n] = 15 - HfN_
is[n] = 16 - ITO_
is[n] = 17 Ц TaN_
is[n] = 17 - TiN_
is[n] = 17 - ZrN_




void libmat(int is, double lam, double *Re_eps, double *Im_eps)
{
	double EPS_r, EPS_i;


	switch (is)
	{
	case 1:
	{
			  EPS_r = Cu(f, 0);
			  EPS_i = Cu(f, 1);
	}

		break;
	case 2:
		Ag_(f, &EPS_r, &EPS_i);

		break;
	case 3:

		Au_(f, &EPS_r, &EPS_i);

		break;
	case 4:
		ZnO_(f, &EPS_r, &EPS_i);
		break;
	case 5:
		Si_(f, &EPS_r, &EPS_i);

		break;
	case 6:

		SiO2_(f, &EPS_r, &EPS_i);

		break;
	case 7:
		HfO_(f, &EPS_r, &EPS_i);
		break;
	case 8:
		AZO_(f, &EPS_r, &EPS_i);
		break;
	case 9:
		EPS_R(f, &EPS_r, &EPS_i);
		break;
	case 10:
		AZO_2(f, &EPS_r, &EPS_i);

		break;
	case 11:

		AZO_wt_0_8_(f, &EPS_r, &EPS_i);

		break;
	case 12:
		AZO_wt_2_0_(f, &EPS_r, &EPS_i);
		break;
	case 13:
		AZO_wt_3_2_(f, &EPS_r, &EPS_i);

		break;
	case 14:

		GZO_(f, &EPS_r, &EPS_i);

		break;
	case 15:
		HfN_(f, &EPS_r, &EPS_i);
		break;
	case 16:
		ITO_(f, &EPS_r, &EPS_i);
		break;
	case 17:
		TaN_(f, &EPS_r, &EPS_i);
		break;
	case 18:
		TiN_(f, &EPS_r, &EPS_i);
		break;
	case 19:
		ZrN_(f, &EPS_r, &EPS_i);
		break;
	}//  end switch 

	*Re_eps = EPS_r;
	*Im_eps = EPS_i;

}     libmat  */




