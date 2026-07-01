/*
Источник: приложение к статье
"Электродинамический анализ оптических поглощающих метаповерхностей,
содержащих диэлектрические слои с плазмонными наночастицами"
(Лерер А.М., Иванова И.Н., Кравченко В.И., Рейзенкинд Я.А.)

Оригинальный файл "COMPOSIT — копия.c" был получен в кодировке Windows-1251
(mojibake при чтении как UTF-8). Ниже — тот же код, перекодированный в UTF-8
без изменения логики. Оригинал лежит в Obsidian-хранилище рядом со статьёй:
PhD/Scientific_study/Scientific_study/paper_analysis/from_Lerer/от Лерера/attachments/COMPOSIT — копия.c

libmat(), divide_(), mult_() определены не в этом файле (объявления и тела
не найдены в доступных исходниках) — они реализуют комплексное деление,
умножение и выбор материала из базы данных диэл. проницаемости соответственно.
Их поведение восстановлено в Python-порте scripts/composite_ema.py по прямому
смыслу вызовов и по формуле Клаузиуса-Моссотти / Максвелла-Гарнетта.
*/

#include "VIBRr.h"

/*
Диэлектрическая проницаемость композита определяетс формулой Клаузиуса-Моссотти.
void COMPOSITE - один сорт н/частиц
void COMPOSITE_3 - 2 сорта н/частиц
*/

void COMPOSITE(double R_eps_m,	double I_eps_m, //Комплексная диэл. проницаемость материала диэлектрика слоя (матрицы)
	double C,//концентрация н/частиц
	int n_mat, //номер из базы данных диэл. проницаемость материала н/частицю 1- Cu, 2-Ag,4-Au и т. д.
	double *R_eps_c, double *I_eps_c)
{
	double Re_eps, Im_eps, rh,ih,rp, ip, u, v;

	libmat(n_mat, f, &Re_eps, &Im_eps);

	u = Re_eps - R_eps_m;
	v = Im_eps - I_eps_m;
	rh = Re_eps + 2.0 * R_eps_m;
	ih = Im_eps + 2.0 * I_eps_m;
	divide_(u, v, rh, ih, &rp, &ip);

	u = 1.0 + 2.0 * C * rp;
	v = 2.0 * C * ip;
	rh = 1.0 - C * rp;
	ih = - C * ip;
	divide_(u, v, rh, ih, &rp, &ip);

	mult_(R_eps_m, I_eps_m, rp,ip, &u, &v);
	*R_eps_c = u;
	*I_eps_c = v;

}//COMPOSITE


void COMPOSITE_3(double R_eps_m, 	double I_eps_m, //Комплексная диэл. проницаемость материала диэлектрика слоя (матрицы)
	double C1, double C2,//концентрации н/частиц
	int n_mat_1, int n_mat_2, //номера из базы данных диэл. проницаемость материала н/частиц
	double* R_eps_c, double *I_eps_c)//Комплексная диэл. проницаемость композита
{
	double Re_eps, Im_eps,  C_[3],rh, ih, u, v, r_CP, i_CP;
	int i,n,I=2;

	C_[1] = C1;
	C_[2] = C2;

	for (r_CP = i_CP = 0.0, i = 1; i <= 2; i++)
	{
		n = (i == 1) ? n_mat_1 : n_mat_2;
		libmat(n, f, &Re_eps, &Im_eps);

		u = Re_eps - R_eps_m;
		v = Im_eps - I_eps_m;
		rh = Re_eps + 2.0 * R_eps_m;
		ih = Im_eps + 2.0 * I_eps_m;
		divide_(u, v, rh, ih, &u, &v);
		r_CP += C_[i] * u;
		i_CP += C_[i] * v;
	}

	u = 1.0 + 2.0 * r_CP;
	v = 2.0 * i_CP;
	rh = 1.0 - r_CP;
	ih = -i_CP;
	divide_(u, v, rh, ih, &u, &v);

	mult_(R_eps_m, I_eps_m, u, v, &u, &v);
	*R_eps_c = u;
	*I_eps_c = v;
}//COMPOSITE 3
