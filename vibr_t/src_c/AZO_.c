#include "VIBRr.h"
void AZO_(double lam_,double *Re_eps,double *Im_eps)
{
	double Lam[20], n[20], k[20], lam;
	int m;
	
	lam = lam_ / 1000.0;

	 Lam[0] = 1.50314;	n[0] = 1.55127;
	 Lam[1] = 1.55924;	n[1] = 1.33414;
	 Lam[2] = 1.63135;	n[2] = 1.03257;
	 Lam[3] = 1.71815;	n[3] = 0.658625;
	 Lam[4] = 1.80227;	n[4] = 0.28468;
	 Lam[5] = 1.87036;	n[5] = -0.053076;
	 Lam[6] = 1.93445;	n[6] = -0.354644;
	 Lam[7] = 1.97317;	n[7] = -0.55971;
	 Lam[8] = 2.02924;	n[8] = -0.837153;
	 Lam[9] = 2.09199;	n[9] = -1.17491;
	 Lam[10] = 2.14405;	n[10] = -1.45235;
	 Lam[11] = 2.20545;	n[11] = -1.81423;
	 Lam[12] = 2.24816;	n[12] = -2.07961;
	 Lam[13] = 2.30957;	n[13] = -2.41737;
	 Lam[14] = 2.38431;	n[14] = -2.88782;
	 Lam[15] = 2.47374;	n[15] = -3.43064;
	 Lam[16] = 2.58585;	n[16] = -4.13028;
	 Lam[17] = 2.6499;	    n[17] = -4.56454;
	 Lam[18] = 2.69527;	n[18] = -4.87817;

	 if (lam < Lam[0] || lam > Lam[18])
	 {
		 printf("\n Re_eps not defined");
		 fprintf(out, "\n Re_eps not defined");
	 
		 fflush(out);
		 exit(1);
	}

	 for (m = 0; m < 18; m++)
	 {
		 if (lam > Lam[m] - 1.0e-11 && lam <= Lam[m + 1])
			 break;
	 }
	 
	 *Re_eps = (n[m + 1] * (lam - Lam[m]) + n[m] * (Lam[m + 1] - lam)) / (Lam[m + 1] - Lam[m]);
	 


	 Lam[0] = 1.5009;	    k[0] = 0.271084;
	 Lam[1] = 1.545;	        k[1] = 0.295181;
	 Lam[2] = 1.65593;     k[2] = 0.361446;
	 Lam[3] = 1.79759;	    k[3] = 0.454819;
	 Lam[4] = 1.93657;	    k[4] = 0.569277;
	 Lam[5] = 2.07554;	    k[5] = 0.698795;
	 Lam[6] = 2.24656;	    k[6] = 0.879518;
	 Lam[7] = 2.3441;	    k[7] = 0.993976;
	 Lam[8] = 2.47368;	    k[8] = 1.16867;
	 Lam[9] = 2.56987;	    k[9] = 1.30422;
	 Lam[10] = 2.70077;	k[10] = 1.50602;

	 if (lam < Lam[0] || lam > Lam[10])
	 {
		 printf("\n Im_eps not defined");
		 fprintf(out, "\n Im_eps not defined");

		 fflush(out);
		 exit(1);
	 }

	 for (m = 0; m < 10; m++)
	 {
		 if (lam > Lam[m] - 1.0e-11 && lam <= Lam[m + 1])
			 break;
	 }

	 *Im_eps = - (k[m + 1] * (lam - Lam[m]) + k[m] * (Lam[m + 1] - lam)) / (Lam[m + 1] - Lam[m]);

	// printf("\n    lam=%6.3lf   *Re_eps=%6.3lf   *Im_eps=%6.3lf ", lam, *Re_eps, *Im_eps);
	// fprintf(out, "\n    lam=%6.3lf   *Re_eps=%6.3lf   *Im_eps=%6.3lf ", lam, *Re_eps, *Im_eps);

 }