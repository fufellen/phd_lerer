#include "VIBRr.h"





void VO2_hot(double x, double *eps1, double *eps2)
{//funcia, kot. vozvraschaet interpolirov-e znacenie eps
	int q=21; // kolichestvo ishodnih tocek
	double  dx = 250.0, e1, e2;
	int i;
	
	double xi[1321]; //massiv dlin voln
	double fs[21][2];//massiv compleksnyh n (n=n'-i*n'')
		
	xi[0]=1000.0;
	xi[1]=1250.0;
	xi[2]=1500.0;
	xi[3]=1750.0;
	xi[4]=2000.0;
	xi[5]=2250.0;
	xi[6]=2500.0;
	xi[7]=2750.0;
	xi[8]=3000.0;
	xi[9]=3250.0;
	xi[10]=3500.0;
	xi[11]=3750.0;
	xi[12]=4000.0;
	xi[13]=4250.0;
	xi[14]=4500.0;
	xi[15]=4750.0;
	xi[16]=5000.0;
	xi[17]=5250.0;
	xi[18]=5500.0;
	xi[19]=5750.0;
	xi[20]=6000.0;



	//n(n'-i*n'') (VO2)
	fs[0][0]=2.2;
	fs[1][0]=2.7;
	fs[2][0]=3.1;
	fs[3][0]=3.65;
	fs[4][0]=4.0;
	fs[5][0]=4.4;
	fs[6][0]=4.95;
	fs[7][0]=5.3;
	fs[8][0] = 5.7;
	fs[9][0] = 6.0;
	fs[10][0] = 6.2;
	fs[11][0] = 6.5;
	fs[12][0] = 6.7;
	fs[13][0] = 7.03;
	fs[14][0] = 7.25;
	fs[15][0] = 7.62;
	fs[16][0] = 7.74;
	fs[17][0] = 7.96;
	fs[18][0] = 8.18;
	fs[19][0] = 8.36;
	fs[20][0] = 8.62;
	

	fs[0][1]=3.5;
	fs[1][1]=4.3;
	fs[2][1]=4.9;
	fs[3][1]=5.4;
	fs[4][1]=5.7;	
	fs[5][1]=6.1;
	fs[6][1]=6.4;
	fs[7][1]=6.7;
	fs[8][1]=7.1;
	fs[9][1]=7.3;
	fs[10][1]=7.7;
	fs[11][1]=8.0;
	fs[12][1]=8.2;	
	fs[13][1] = 8.3;
	fs[14][1] = 8.65;
	fs[15][1] = 8.8;
	fs[16][1] = 9.07;
	fs[17][1] = 9.24;
	fs[18][1] = 9.46;
	fs[19][1] = 9.7;
	fs[20][1] = 9.85;
	for (i = 0; i < q; i++)
	{
		if (fabs(x - xi[i]) < 1.0e-7)
		{
			e1 = fs[i][0];
			e2 = fs[i][1];
			break;
		}
		else
		{
			if (x > xi[i] && x < xi[i+1])
			{
				e1 = (fs[i + 1][0] * (x - xi[i]) - fs[i][0] * (x - xi[i + 1])) / dx;
				e2 = (fs[i + 1][1] * (x - xi[i]) - fs[i][1] * (x - xi[i + 1])) / dx;
				break;
			}
		}
	}//i

	//printf("\n  l=%6.3lf   n_1=%6.3lf    n_2=%6.3lf", x, e1, e2);
	*eps1 = e1 * e1 - e2 * e2;
	*eps2 = -2.0 * e1 * e2;
	

}

