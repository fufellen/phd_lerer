
#include "VIBRr.h"

void VVOD ( void )
{
	int n,m;
	double u;
	// Решетка диэлектрических эллиптических цилиндров
	// или эллиптических отверстий в диэлектрике.  10.11.10
	//Косая решетка из эллиптических цилиндров. База данных для металлов и диэлектриков. Композитные материалы



		
	/*
	Параметры цикла по длине волны в произвольных единицах
	*/

	fSkipComment( in );
	fscanf (in,"%lf %lf %lf",&AF[1],&AF[2],&AF[3]);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	AF[1] += 0.0001;

	// периоды структуры и тангнс угла наклона решетки
	fSkipComment( in );
	fscanf (in,"%lf %lf %lf",&dx,&dy,&tg_pci);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!2


	
	// Obshee chislo sloev - Nsl + 1. Numeratcia sverhu: 
	//pervii polubeskonechnii - otkuda padaet volna (vozduh), Nsl + 1 - polubeskonechnaia podlojka.
	
	/*
	Нумерация диэлектриков
	0 - неоднородность.n_[0], k_[0] - действительная и мнимая части показателя преломления внутренней части неоднородности
	n_in, k_in - действительная и мнимая части показателя преломления внешней части неоднородности

	1 - диэлектрик над решеткой.Из него падает волна.Этот диэлектрик может быть :
	полубесконечный, для этого нужно положить H_[0] = 0;
	б) если положить H_[0] = 1, то диэлектрик двухслойный, верхний слой - полубесконечный воздух,
	затем слой с показателем преломления n_[1], коэффициенты преломления и отражения рассчитываются
	по формуле Френеля.В этом случае толщина слоя не влияет на характеристики.

	2 - слой, в котором находится неоднородность.Его толщина H_[1].
	3, 4, 5 - подложка.Толщина 3 - го слоя - H_[2], 4 - г - H_[3].Слой 5 - полубесконечный.

	Если is[n] = 0, то считывается из файла с ч.д
	{
	n_[n], k_[n] - действительная и мнимая части показателя преломления.
	k_[1] = 0 ОБЯЗАТЕЛЬНО!!!!!!
	}
	Если is[n] != 0, то определяются в базе данных.
	is[n] = 1 медь
	is[n] = 2 - серебро
	is[n] = 3 - золото
	is[n] = 4 - ZnO
	is[n] = 5 - Si
	is[n] = 6 - SiO2
	is[n] = 7 - HfO
	is[n] = 8 - AZO

	см библиотеку материалов!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	*/
	/*switch (is_l[n][m])
	{
	case 1:
	{
	EPS_l_r[n][m] = Cu(f, 0);
	EPS_l_i[n][m] = Cu(f, 1);
	}

	break;
	case 2:
	Ag_(f, &EPS_l_r[n][m], &EPS_l_i[n][m]);

	break;
	case 3:

	Au_(f, &EPS_l_r[n][m], &EPS_l_i[n][m]);

	break;
	case 4:
	ZnO_(f, &EPS_l_r[n][m], &EPS_l_i[n][m]);
	break;
	case 5:
	Si_(f, &EPS_l_r[n][m], &EPS_l_i[n][m]);
	break;
	case 6:
	SiO2_(f, &EPS_l_r[n][m], &EPS_l_i[n][m]);
	break;
	case 7:
	HfO_(f, &EPS_l_r[n][m], &EPS_l_i[n][m]);
	break;
	case 8 и 10 - более полная:
	AZO_(f, &EPS_l_r[n][m], &EPS_l_i[n][m]);
	break;
	}  end switch */



	// N_hol количество слоев, в которых расположена неоднородность. Нумерация этих слоев 2, .... N_hol .
	//N_ell - количество эллипсов
	fSkipComment( in );
	fscanf (in,"%d %d %d",&Nsl,&N_hol,&N_ell);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	AM1( );//new
	AM2( );//new	
	for ( N_La=0, n=1; n <= N_hol; n++ )
	{
		/*
		внешние размеры неоднородности
		/A_t[n], B_t[n] - размеры осей эллипса на верхней границе n-го отверстия (нумерация отверстий сверху вниз)
		A_b[n], B_b[n] - размеры осей00ллипса на нижней  границе n-го отверстия 
		N_l[n] - количество лайнеров в неоднородности в слое с номером n
		fi_0 - угол поворота оси A к оси X
		*/
		fSkipComment( in );
		//new
		fscanf (in,"%lf %lf %lf %lf %d %lf",&A_t[n],&B_t[n],&A_b[n],&B_b[n],&N_l[n],&fi_0[n]);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		if (N_La < N_l[n])
			N_La = N_l[n];
		//new end
	}

	AM3( );//new

	//параметры лайнеров
	//  N_l[n] - количество лайнеров в неоднородности в слое с номером n.  N_l[n]<6
	for (bum = 0, u = 0, n = 1; n <= N_hol; n++)
	{
		for ( m=1;m <= N_l[n]; m++ )
		{
			/*
			параметры лайнеров
			A_t_l[n][m], B_t_l[n][m] - оси m-го лайнера на верхней границе n-го отверстия
			A_b_l[n][m], B_b_l[n][m] - оси m-го лайнера на нижней  границе n-го отверстия 
			
			n_l[n][m], k_l[n][m] - действительная и мнимая часть комплексного показателя преломления
			m-го лайнера n-го отверстия или основы композита
			
			В ПРОГРАММЕ ПОЛАГАТСЯ, что внутри слоя 
			A_t_l[n][m] / B_t_l[n][m] != A_t[n] / B_t[n],  A_b_l[n][m] / B_b_l[n][m] != A_b[n] / B_b[n]
			
			is_l[n][m] номер материала слоя в неоднородности
			comp_l[n][m] - номер материала композита слоя в неоднородности
			CONCEN_l - концентрация наночастиц  в неоднородности
			*/

			fSkipComment(in);
			fscanf(in, "%d %d %lf", &is_l[n][m], &comp_l[n][m], &CONCEN_l[n][m]);
			//fscanf(in, "%d", &is_l[n][m]);
			bum += is_l[n][m];
			
			if (m == N_l[n])
			{
				A_t_l[n][m] = A_t[n];
				A_b_l[n][m] = A_b[n];
				
				B_t_l[n][m] = B_t[n];
				B_b_l[n][m] = B_b[n];
				
				if (!is_l[n][m])
				{
					fSkipComment( in );
					fscanf (in,"%lf %lf",&n_l[n][m], &k_l[n][m]);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				}
			}//m == N_l[n]
			else
			{
				if (!is_l[n][m])
				{
					//fSkipComment( in );
					fscanf(in, "%lf %lf %lf %lf %lf %lf", &A_t_l[n][m], &B_t_l[n][m], &A_b_l[n][m], &B_b_l[n][m], &n_l[n][m], &k_l[n][m]);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				}
				else
				{
					//fSkipComment( in );
					fscanf(in, "%lf %lf %lf %lf", &A_t_l[n][m], &B_t_l[n][m], &A_b_l[n][m], &B_b_l[n][m]);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				}
			}//!m == N_l[n]

			u += (A_t_l[n][m] / A_t[n] - B_t_l[n][m] / B_t[n]) + ( A_b_l[n][m] / A_b[n]- B_b_l[n][m] / B_b[n]);
	
			if (!is_l[n][m])
			{
				EPS_l_r_[n][m] = n_l[n][m] * n_l[n][m] - k_l[n][m] * k_l[n][m];
				EPS_l_i_[n][m] = -2.0 * n_l[n][m] * k_l[n][m];
			}
		}//m
	}//n

	Ekv = (fabs(u) > 0.0001) ? 0 : 1;
	// координаты центров эллипсов
	for ( n=1; n <= N_ell; n++ )
	{
		fSkipComment( in );
		fscanf (in,"%lf %lf",&XC[n],&YC[n]);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	}
	
	/*
		Нумерация диэлектриков
		
		1 - диэлектрик над решеткой. Из него падает волна. Этот диэлектрик может быть:
				полубесконечный, для этого нужно положить H_[0] = 0;
				
		2 - N_hol +1  слои с неоднородностями. И толщина H_[1].
		Общее число слоев Nsl+1
		Если is[n] = 0, то считывается из файла с ч.д
		{
			n_[n], k_[n] - действительная и мнимая части показателя преломления.
			k_[1] = 0 ОБЯЗАТЕЛЬНО!!!!!!
		}
		
		comp_[n] - номер материала композита слоя
		CONCEN - концентрация наночастиц
	 */

	H_s = 0.0;
	for (bum1 = 0, n = 1; n <= Nsl + 1; n++)
	{ // loop over 6 materials
		fSkipComment( in );
		fscanf(in, "%d %d %lf", &is[n], &comp_[n], &CONCEN[n]);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		bum1 += is[n];
		
		//if (!is[n])
		if (n <=1 || n == Nsl+1)
		{
			//fSkipComment( in );
			if (!is[n])
				fscanf (in,"%lf %lf ",&n_[n],&k_[n]);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
		}
		else
		{
			//fSkipComment( in );
			if (!is[n])
				fscanf (in,"%lf %lf %lf ",&H_[n-1],&n_[n],&k_[n]);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			else
				fscanf(in, "%lf", &H_[n - 1]);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
			H_s += H_[n-1];
		}
		
		bu_p[n] = (k_[n] > 1.0e-10) ? 1 : 0;
		
		if (!is[n])
		{
			EPS_r[n] = EPS_r_[n] = n_[n] * n_[n] - k_[n] * k_[n];//!!!!!!!!!!!!
			EPS_i[n] = EPS_i_[n] = -2.0 * n_[n] * k_[n];
		}

	} // // loop over 6 materials
		
	// 2*t_ - общая толщина слоев с неоднородностью

	for ( t_ = 0.0, n=1; n <= N_hol; n++ )
	{
		h_[n] = 0.5 * H_[n];
		t_ += h_[n];
	}

		
	/*
	 Углы падения и расположения неоднородностей в системе координат, у которой 
		ось Z перпендикулярна подложке, 
		ось X Y - в плоскости слоев
		
		tet fi - углы падения внешней волны сферической системы координат в градусах
	*/

	/*внешнее поле exp(i(kx*x+ky*y+kz*z)) пока из первой области.
	
	  Угол tet - отсчитывается от нормали к подлжке -от оси Z.
	Угол fi - отсчитывается от оси X.
	*/

	
	fSkipComment( in );
	fscanf (in,"%lf  %lf",&TET,&FI);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	TET_r = TET * Pi / 180.0;
	FI_r = FI * Pi / 180.0;

	//FI_r_ = FI_r;
	FI_r_ = Pi * 0.5;
	//FI_r_ = 0.0;
	/*
	SIM_r=0 Частичный учет симметрии матрицы - только по M_z. Кроме того, принудительно задается M_r = M_fi
	SIM_r=-1 Частичный учет симметрии матрицы - только по M_z. Цикл по M_r работает
	*/
	/*
	M_z,M_fi,M_r - число  базисных функций по z, fi, r
		*/

	SIM_r = - 1;
	
	//M_z coefficient for numkber of Z-axis base functions;
	// AM_z loop start[1], step[2], finish[3] for M_r optimization
	fSkipComment( in );
	fscanf (in,"%d %d %d",&AM_z[1],&AM_z[2],&AM_z[3]);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	//M_fi number of angular base functions;
	// AM_fi loop start[1], step[2], finish[3] for M_fi optimization
	fSkipComment( in );
	fscanf (in,"%d %d %d",&AM_fi[1],&AM_fi[2],&AM_fi[3]);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	//M_r number of radial base functions;
	// AM_r loop start[1], step[2], finish[3] for M_r optimization
	fSkipComment( in );
	fscanf (in,"%d %d %d",&AM_r[1],&AM_r[2],&AM_r[3]);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	


	/*
	N_g - число членов в функции Грина.
	*/
	// N_g half truncation number for q and p Flouqeut harmonics
	// A.L. N_g for p and q should be separated!!!
	// AN_g loop start[1], step[2], finish[3] for N_g optimization
	fSkipComment( in );
	fscanf (in,"%d %d %d",&AN_g[1], &AN_g[2], &AN_g[3]);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	// N_az, N_rad - порядки квадратур по азимуту и радиусу для вычисления интегралоа для лайнеров
	// у которых A_t_l[n][m] / B_t_l[n][m] != A_t[n] / B_t[n],  A_b_l[n][m] / B_b_l[n][m] != A_b[n] / B_b[n]
	fSkipComment( in );
	fscanf (in,"%d %d ",&N_az, &N_rad);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	// 	ft_l число << 1, влияющее на точность вычисдения ядра
	fSkipComment( in );
	fscanf (in,"%lf",&ft_l);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	//Матричные элементы - сумма членов двойных рядов, ограниченны по максимльному радиусу по ро - Rmax. 
	//Этот радиус определяется значениямт Ng. 
	//Ряды для однородной части функции Грина сходятся гораздо лучше, чем рядв для однородной части ФГ. Поэтому эти ряды суммируются в круге радиуса
	//Rmax / Rad_ro.  Rad_ro>=1. Рекомендуемое значение  5>=Rad_ro>=2
	fSkipComment( in );
	fscanf (in,"%lf",&Rad_ro);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	

//Конец ввода числовых данных
	u = 0.001;
	BU_SIM = (fabs(FI) < u || fabs(FI - 90.0) < u || fabs(FI + 90.0) < u || fabs(FI - 180.0) < u || fabs(FI + 180.0) < u || fabs(FI - 270.0) < u) ? 1 : 0;
		
//____________________________________________________________________________________________________________	

	printf ( "\n\n    N_az =%3d    N_rad=%3d   ft_l=%6.3e",N_az,N_rad,ft_l);
	fprintf ( out,"\n\n    N_az =%3d    N_rad=%3d   ft_l=%6.3e",N_az,N_rad,ft_l);

	printf ( "\n\n    TET =%6.3lf    FI=%6.3lf",TET,FI);
	fprintf ( out,"\n\n    TET =%6.3lf    FI=%6.3lf",TET,FI);	
	
	printf ( "\n\n  dx =%6.3lf  dy =%6.3lf",dx,dy);
	fprintf ( out,"\n\n  dx =%6.3lf  dy =%6.3lf",dx,dy);
	
	for ( n=1; n <= N_hol; n++ )
	{
		//new
		printf ( "\n\n  n=%2d Ax_t[n]=%6.3lf Ay_t[n]=%6.3lf Ax_b[n]=%6.3lf Ay_b[n]=%6.3lf fi_0[n]=%5.2lf\n",n,A_t[n],B_t[n],A_b[n],B_b[n], fi_0[n]);
		fprintf ( out,"\n\n  n=%2d Ax_t[n]=%6.3lf Ay_t[n]=%6.3lf Ax_b[n]=%6.3lf Ay_b[n]=%6.3lf fi_0[n]=%5.2lf\n",n,A_t[n],B_t[n],A_b[n],B_b[n], fi_0[n]);

		fi_0[n] *= Pi /180.0;

		c_fi0[n] = cos(fi_0[n]);
		s_fi0[n] = sin(fi_0[n]);
		//new end
	}
	for ( n=1; n <= Nsl+1; n++ )
	{
		if (n > 1 && n < Nsl+1)
		{
			if (is[n])
			{
				printf ( "\n  n=%1d metall is[n]=%1d  H[%1d] =%6.3lf",n,is[n],n,H_[n-1]);
				fprintf ( out,"\n  n=%1d metall is[n]=%1d  H[%1d] =%6.3lf",n,is[n],n,H_[n-1]);
			}
			else
			{
				printf ( "\n  n[%1d]=%8.5lf  k[%1d] =%8.5lf  H[%1d] =%6.3lf",n,n_[n],n,k_[n],n-1,H_[n-1]);
				fprintf ( out,"\n  n[%1d]=%8.5lf  k[%1d] =%8.5lf  H[%1d] =%6.3lf",n,n_[n],n,k_[n],n-1,H_[n-1]);
			}
			if (n <= 1 + N_hol)
			{
				printf ( "  Hole");
				fprintf ( out,"   Hole");
			}
		}
		else
		{
			if (is[n])
			{
				printf ( "\n n=%1d metall is[n]=%1d",n,is[n]);
				fprintf ( out,"\n n=%1d metall is[n]=%1d",n,is[n]);
			}
			else
			{
				printf ( "\n  n[%1d]=%8.5lf  k[%1d] =%8.5lf",n,n_[n],n,k_[n]);
				fprintf ( out,"\n  n[%1d]=%8.5lf  k[%1d] =%8.5lf",n,n_[n],n,k_[n]);
			}
		}
	}//n

	for ( n=1; n <= N_hol; n++ )
		for ( m=1;m <= N_l[n]; m++ )
		{
			if (!is_l[n][m])
			{
				printf ( "\n\n  n=%2d  m=%2d \n   A_t[n][m]=%5.3lf  B_t[n][m]=%5.3lf  A_b[n][m]=%5.3lf  B_b[n][m]=%5.3lf\n   n_l[n][m]=%5.3lf  k_l[n][m]=%5.3lf",
					n,m,A_t_l[n][m],B_t_l[n][m], A_b_l[n][m],B_b_l[n][m], n_l[n][m], k_l[n][m]);
				fprintf ( out,"\n\n  n=%2d  m=%2d\n A_t[n][m]=%5.3lf B_t[n][m]=%5.3lf A_b[n][m]=%5.3lf B_b[n][m]=%5.3lf\n   n_l[n][m]=%5.3lf  k_l[n][m]=%5.3lf",
					n,m,A_t_l[n][m],B_t_l[n][m], A_b_l[n][m],B_b_l[n][m], n_l[n][m], k_l[n][m]);
			}
			else
			{
				printf("\n\n  n=%2d  m=%2d \n   A_t[n][m]=%5.3lf  B_t[n][m]=%5.3lf  A_b[n][m]=%5.3lf  B_b[n][m]=%5.3lf   metall is[n][m]=%1d",
					n, m, A_t_l[n][m], B_t_l[n][m], A_b_l[n][m], B_b_l[n][m], is_l[n][m]);
				fprintf(out, "\n\n  n=%2d  m=%2d\n A_t[n][m]=%5.3lf B_t[n][m]=%5.3lf A_b[n][m]=%5.3lf B_b[n][m]=%5.3lf   metall is[n][m]=%1d",
					n, m, A_t_l[n][m], B_t_l[n][m], A_b_l[n][m], B_b_l[n][m], is_l[n][m]);
			}
		}

	printf ( "\n");
	fprintf ( out,"\n");
	for ( n=1; n <= N_ell; n++ )
	{
		printf ( "\n  XC[%1d]=%5.3lf  YC[%1d]=%5.3lf",n,XC[n],n,YC[n]);
		fprintf ( out,"\n XC[%1d]=%5.3lf  YC[%1d]=%5.3lf",n,XC[n],n,YC[n]);
	}
		
}   /*  end VVOD  */


	  /*
		A, B -оси внешненго эллипса, образующего неоднородность. a_, b_ -его полуоси
	

	*/

	/*
		 Углы падения и расположения вибраторов в системе координат, у которой 
			ось Z перпендикулярна подложке, вдоль вибратора
			ось X перпендикулярна вибратору.
			ось Y перпендикулярна вибратору.
			
		*/
	

	
