
#include "VIBRr.h"

void F_df_2(int sp, //0-s-pol, 1-p-pol
	double *eps, double *mu, double *R, double *T)
{
	double   gam, kk, x, y, x1,y1,Rek, Imk, kh[3], Rec, Imc, Res, Ims, 	del[3],u, v, Imkap, Rekap,gam2,rt,it,
		Imkap_mu, Rekap_mu, Rek_mu, Imk_mu;

	if (kbhit())
		if (getch() == ESCape)
		{
			fflush(out);
			exit(1);
		}
	
	
	kk = kx * kx;
	gam2 = k2 - kk;
	gam = sqrt(gam2);

	mult(eps[1], eps[2], mu[1], mu[2], &x, &y);
	Rekap =  k2 * x - kk;//kappa^2
	Imkap = k2 * y;

	if (h <= 10000.0)
	{
		if (!sp)
		{
			divide(Rekap, Imkap, mu[1], mu[2], &Rekap_mu, &Imkap_mu);
			divide(Rekap_mu, Imkap_mu, mu[1], mu[2], &Rekap_mu, &Imkap_mu);// kappa^2 s chertoy
		}
		else
		{
			divide(Rekap, Imkap, eps[1], eps[2], &Rekap_mu, &Imkap_mu);
			divide(Rekap_mu, Imkap_mu, eps[1], eps[2], &Rekap_mu, &Imkap_mu);// kappa^2 s chertoy
		}

		compl_sqrt(Rekap_mu, Imkap_mu, &Rek_mu, &Imk_mu);

		compl_sqrt(Rekap, Imkap, &Rek, &Imk);
		kh[1] = Rek * h;
		kh[2] = Imk * h;

		compl_cos_sin(kh[1], kh[2], &Rec, &Imc, &Res, &Ims);

		mult(Rec, Imc, Rek_mu, Imk_mu, &x, &y);// T
		mult(x, y, 0.0, 2.0 * gam, &rt, &it);

		x = Rekap_mu + gam2;
		y = Imkap_mu;
		mult(x, y, Res, Ims, &x1, &y1);
		del[1] = rt - x1;
		del[2] = it - y1;

		mult(Rek_mu, Imk_mu, 0.0, 2.0 * gam, &rt, &it);
		divide(rt, it, del[1], del[2], &u, &v);
		*T = (u * u + v * v);

		x = Rekap_mu - gam2;
		y = Imkap_mu;
		mult(x, y, Res, Ims, &x, &y);
		divide(x, y, del[1], del[2], &u, &v);

		*R = (u * u + v * v);
	}//h <= 10000.0
	else
	{
		compl_sqrt(Rekap, Imkap, &Rek, &Imk);
		rt = Rek / gam;
		if (sp)
		{
			divide(Rek, Imk, eps[1], eps[2], &x, &y);
			rt /= eps[1];
		}
		else
		{
			divide(Rek, Imk, mu[1], mu[2], &x, &y);
			rt/= mu[1];
		}
		x /= gam;
		y /= gam;

		divide(1.0 - x, -y, 1.0 + x, y,  &u, &v);
		*R = (u * u + v * v);
		divide(2.0, 0.0, 1.0 + x, y, &u, &v);

		*T = (fabs(EPS[2]) < 0.000000001 && fabs(mu[2]) < 0.000000001) ? (u * u + v * v) * rt: 0.0;
	}//h > 10000.0
	
}//F_df
