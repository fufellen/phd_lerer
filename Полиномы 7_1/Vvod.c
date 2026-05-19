
#include "VIBRr.h"

void VVOD ( void)
{

	int m;
	// Нахождение ДП и tg delta диэлектрической пластины с потерями.  4.12.25
	
	//TEST = 1 отражение от слоя в вакууме
	//TEST = 2 отражение от слоя на подложке
	fSkipComment(in);
	fscanf(in, "%d ", &TEST);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	
	//TET	 угол падения в градусах, 
	
	fSkipComment( in );
	fscanf (in,"%lf ",&TET);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	//h,  l - толщина и ширина пластины, d - период
	// Если частота в ГГц, то в размеры в мм.
	// Если частота в ТГц, то в размеры в мкм.
	fSkipComment(in);
	fscanf(in, "%lf %lf %lf", &d, &h, &l);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	TET *= Pi / 180.0;
	//AF[1:3] частоты
	fSkipComment(in);
	fscanf(in, "%lf %lf  %lf", &AF[1], &AF[2],&AF[3]);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	//EPS[1:2],MU[1:2] - неоднородности
	
	fSkipComment(in);
	fscanf(in, "%lf %lf", &EPS[1], &EPS[2]);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	EPS[2] = -EPS[2];

	fSkipComment(in);
	fscanf(in, "%lf %lf", &MU[1], &MU[2]);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	MU[2] = -MU[2];
	
	//eps[2:3], mu__[2:3] - слой толщина h2 и подложка
	fSkipComment(in);
	fscanf(in, "%lf ", &h2);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	
	for (m = 2; m <= 3; m++)
	{
		fSkipComment(in);
		fscanf(in, "%lf %lf", &eps[m], &mu__[m]);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	}
		
	//M_X,M_Y - число  базисных функций по x, y
	fSkipComment(in);
	fscanf(in, "%d %d %d", &AM_X[1], &AM_X[2], &AM_X[3]);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	fSkipComment(in);
	fscanf(in, "%d %d %d", &AM_Y[1], &AM_Y[2], &AM_Y[3]);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	/*
	N_g - число членов в функции Грина.
	*/
	fSkipComment(in);
	fscanf(in, "%d %d %d", &AN_g[1], &AN_g[2], &AN_g[3]);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	//PRINT = 0 - не печатать в файл
	//PRINT = 1 - печать в файл
	fSkipComment(in);
	fscanf(in, "%d ", &PRINT);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	printf("\n  d=%6.4lf     l=%6.4lf    h=%6.4lf \n eps: %6.4lf   i*%6.4lf\n mu: %6.4lf   i*%6.4lf\n", d, l, h, EPS[1], EPS[2], MU[1], MU[2]);
	if (PRINT)
	fprintf(out, "\n  d = %6.4lf     l = %6.4lf    h = %6.4lf \n eps : %6.4lf   i*%6.4lf\n mu : %6.4lf   i*%6.4lf\n", d, l, h, EPS[1], EPS[2], MU[1], MU[2]);
	SX = (AM_Y[2] < 0) ? 1 : 0;
	//printf("\n mu:  %6.4lf  i*%6.4lf", MU[1], MU[2]);

	for (m = 2; m <= 3; m++)
	printf("\n  m = %1d eps=%6.4lf     mu=%6.4lf", m,eps[m],mu__[m]);

	printf("\n  ");

}   /*  end VVOD  */


	  
