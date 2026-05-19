
#include "VIBRr.h"


int main ( void )
{
	int var, Nvar;
	double  R,T,P;
		
	//  К определению расстояния.  30.03.21

	printf("\n Enter the quantity of variants - ");
	scanf("%d",&Nvar);
	//Nvar = 1;

	printf("\n Enter the name for input file - ");

	scanf("%s",inf);


	if ((in = fopen(inf,"r")) == NULL)
	{
		printf("\n Error! Input file was not found! \n");
		exit(1);
	}
	printf("\n Enter the name of the file for results - ");

	scanf("%s",outf);

	if ((out = fopen(outf,"a")) == NULL)
	{
		printf("\n The target file is not opened! \n");
		exit(1);
	}

	/*
	printf("\n ");

	
	printf("\n Will you plot graphs F(frequency)?(Y/N) - ");
	graf2 = yesno();

	if (graf2)
	{
		printf("\n Enter a name of a file for record of the table for construction of diagrams - ");
		scanf("%s",out2f);
		if ((out2 = fopen(out2f,"a")) == NULL)
		{
			printf("\n  The target file is not opened! \n");
			exit(1);
		}
		
				
	}//graf2*/

	Pi = 4.0 * atan(1.0);
	Pi_2 = 2.0 * Pi;
	Pi05 = Pi * 0.5;
	G_R = Pi / 180.0;

	var = 0;
	start = clock();
	
	
	
/*--------------------------------------------------------------------------*/
	do
	{	
		str = clock();
		var++;


		VVOD(); // input 
		//test
		/*{
			double z=3.0, I, K, mn;

			Imod_K(3, z, &I, &K);
			
			K = I_ex(3,  z) * exp(z);
			printf("\n\n i=%9.7e    k=%9.7e ",  I, K);*/
			
			/*I = Pn(122, z);
			K = Pn(123, z);
			printf("\n\n i=%9.7e    k=%9.7e ", I,  K);
			
			mn = sqrt(0.5 * Pi);
			
			Imod_K(0, z, &I, &K);
			I = J(0, z);
			K = N(0, z);
			printf("\n\n i=%9.7e    k=%9.7e ", mn *  I, mn * K);
			I = sqrt(2.0 / Pi) / z * sin(z);
			K = -sqrt(2.0 / Pi) / z * cos(z);
			printf("\n\n i=%9.7e    k=%9.7e ", mn *  I, mn * K);

			Imod_K(0, z, &I, &K);
			printf("\n\n i=%9.7e    k=%9.7e ", mn *  I, mn * K);
			I = sqrt(2.0 / Pi) / z * sinh(z);
			K = 1.0 / sqrt(2.0 / Pi) / z * exp(-z);
			printf("\n\n i=%9.7e    k=%9.7e ", mn *  I, mn * K);

			I = J(1, z);
			K = N(1, z);

			printf("\n\n i=%9.7e    k=%9.7e ", mn *  I, mn * K);
			I = J(2, z);
			K = N(2, z);

			printf("\n\n i=%9.7e    k=%9.7e ", mn *  I, mn * K);

			Imod_K(1, z, &I, &K);
			printf("\n\n i=%9.7e    k=%9.7e ", mn *  I, mn * K);

			Imod_K(2, z, &I, &K);
			printf("\n i=%9.7e    k=%9.7e ", mn *  I, mn * K);*/
		//	exit(1);
		//}//test
	
		
		for (f = AF[1]; f <= AF[3]; f += AF[2])
		{
			k = f * Pi / 149.9;
			k2 = k * k;

			kx = k * sin(TET);

			if (TEST == 1)
			{
				F_df_2(1, EPS, MU, &R, &T);// Аналитическое решение для  слоя в вакууме- p
				P = 1.0 - R - T;
				printf("\n p(anal)  f=%6.4lf   R=%6.4e    T=%6.4e    P=%6.4e ", f, R, T, P);
				if (PRINT)
					fprintf(out, "\n p(anal)  f=%6.4lf   R=%6.4e    T=%6.4e    P=%6.4e ", f, R, T, P);

				F_df_2(0, EPS, MU, &R, &T);// Аналитическое решение для слоя в вакууме - s
				P = 1.0 - R - T;
				printf("\n s(anal)  f=%6.4lf   R=%6.4e    T=%6.4e    P=%6.4e ", f, R, T, P);
				if (PRINT)
					fprintf(out, "\n s(anal)  f=%6.4lf   R=%6.4e    T=%6.4e    P=%6.4e ", f, R, T, P);
												
			}///TEST=1

			if (TEST == 2)
			{
				double rR, iR, rTA, iTA;
				
				epsG[2] = EPS[1];
				epsG[3] = eps[3];
				muG[2] = MU[1];
				muG[3] = mu__[3];
				h2G = h2 +  h;

				E_ext(1, 1,&rR, &iR, &rTA, &iTA);// Аналитическое решение для слоя на подложке - p
				R = rR * rR + iR * iR;
				T = rTA * rTA + iTA * iTA;
				P = 1.0 - R - T;
				printf("\n p(anal)  f=%6.4lf   R=%6.4e    T=%6.4e    P=%6.4e ", f, R, T, P);
				if (PRINT)
					fprintf(out, "\n p(anal)  f=%6.4lf   R=%6.4e    T=%6.4e    P=%6.4e ", f, R, T, P);

				E_ext(0, 1,&rR, &iR, &rTA, &iTA);// Аналитическое решение для слоя на подложке - s
				R = rR * rR + iR * iR;
				T = rTA * rTA + iTA * iTA;
				P = 1.0 - R - T;
				printf("\n s(anal)  f=%6.4lf   R=%6.4e    T=%6.4e    P=%6.4e ", f, R, T, P);
				if (PRINT)
					fprintf(out, "\n s(anal)  f=%6.4lf   R=%6.4e    T=%6.4e    P=%6.4e ", f, R, T, P);

				//exit(1);
			}///TEST/ 

			for (N_g = AN_g[1]; N_g <= AN_g[3]; N_g += AN_g[2]) 
			for (M_X = AM_X[1]; M_X <= AM_X[3]; M_X += AM_X[2])
			{
				if (SX)
					AM_Y[1] = AM_Y[2] = AM_Y[3] = M_X;
				for (M_Y = AM_Y[1]; M_Y <= AM_Y[3]; M_Y += AM_Y[2])
				{
					/*F_df_IE_MK_test(EPS, MU, &R, &T);
					fflush(out);
					exit(1);*/
					MS = M_X * M_Y;
					printf("\n\n  N_g=%2d  M_X=%2d     M_Y=%2d     MS=%3d", N_g, M_X, M_Y, MS);
					if (PRINT)
						fprintf(out, "\n\n  N_g=%2d  M_X=%2d     M_Y=%2d     MS=%3d", N_g, M_X, M_Y, MS);
					{
						VAR = 0;
						printf("\n    ");
						if (PRINT)
							fprintf(out, "\n    ");

						PP = 1;
						F_df_IE_MK(1, &R, &T);//2.	Решаем СЛАУ методом коллокации
						P = 1.0 - R - T;
						printf("     p(numer)  R=%6.4e    T=%6.4e    P=%6.4e",  R, T, P);
						if (PRINT)
							fprintf(out, " p(numer)       R=%6.4e    T=%6.4e    P=%6.4e", R, T, P);
						PP = 0;
						F_df_IE_MK(0, &R, &T);//2.	Решаем СЛАУ методом коллокации
						P = 1.0 - R - T;
						printf("        s(numer)      R=%6.4e    T=%6.4e   P=%6.4e",  R, T, P);
						if (PRINT)
							fprintf(out, "      s(numer)   R=%6.4e   T=%6.4e   P=%6.4e", R, T, P);
					}
					
				}//M_Y
			}//M_x

		}//f
		

		
		tim ( str);
		printf ( "\n\n\n__________________________________________________________________");
		fprintf ( out,"\n\n\n__________________________________________________________________");

	} while( var <  Nvar );

	
	
	if ( Nvar > 1 )
		tim ( start);
	return ( 0 );
}   /*  end main  */


