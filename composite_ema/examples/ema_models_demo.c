/*
  Пример использования и автономная проверка C-моделей эффективной среды
  (../source_c/COMPOSIT.c + ../source_c/COMPOSITE_MODELS.c).

  Демо собирает РЕАЛЬНЫЕ функции моделей (включает оба .c ниже) и прогоняет их
  сквозным путём, как в программе Лерера: задаём длину волны в глобальной f,
  выбираем материал частиц по номеру n_mat через libmat(), считаем eps_eff
  тремя моделями:
     COMPOSITE()      — Клаузиус-Моссотти / Максвелл-Гарнетт (формула (1));
     BRUGGEMAN()      — модель Бруггемана (формула (4));
     COMPOSITE_MLWA() — Максвелл-Гарнетт с размерной поправкой MLWA (R_p).

  Материальная база демо — таблицы Джонсона-Кристи 1972 [12] для Au и Cu
  (49 узлов, конвенция eps = eps' - i*eps''), та же, что в Python-порте.
  n_mat: 4 = Au, 1 = Cu (как в комментарии ../source_c/COMPOSIT.c).

  Прогон печатает положение и амплитуду пика |Im eps_eff| и сравнивает с
  таблицей 2 статьи. Сборка/запуск — см. run_examples.bat / run_examples.ps1
  или examples/README.md.

  ВНИМАНИЕ: divide_/mult_/libmat/f и заголовок VIBRr.h здесь — минимальные
  заглушки ТОЛЬКО для автономной сборки демо. В проекте Лерера используются
  его собственные версии.
*/

#include "VIBRr.h"

/* --- глобальная длина волны (нм), которую читают обёртки моделей --- */
double f = 0.0;

/* --- комплексные помощники (как в newLIb/MUL_DIV.C) --- */
void mult_(double x1, double y1, double x2, double y2, double *u, double *v)
{
	*u = x1 * x2 - y1 * y2;
	*v = x1 * y2 + y1 * x2;
}
void divide_(double x, double y, double x1, double y1, double *u, double *v)
{
	double z = 1.0 / (x1 * x1 + y1 * y1);
	*u = ( x * x1 + y * y1) * z;
	*v = (-x * y1 + y * x1) * z;
}

/* --- материальная база: Johnson & Christy 1972, Au и Cu (eps = eps' - i eps'') --- */
static const double AU_XI[49] = {187.9, 191.6, 195.3, 199.3, 203.3, 207.3, 211.9, 216.4, 221.4, 226.2, 231.3, 237.1, 242.6, 249, 255.1, 261.6, 268.9, 276.1, 284.4, 292.4, 300.9, 310.7, 320.4, 331.5, 342.5, 354.2, 367.9, 381.5, 397.4, 413.3, 430.5, 450.9, 471.4, 495.9, 520.9, 548.6, 582.1, 616.8, 659.5, 704.5, 756, 821.1, 892, 984, 1088, 1216, 1393, 1610, 1937};
static const double AU_RE[49] = {0.2271, 0.2952, 0.2925, 0.2039, 0.1382, -0.0104, -0.1325, -0.2338, -0.3463, -0.4155, -0.551, -0.6169, -0.7445, -0.8913, -1.0804, -1.2365, -1.3464, -1.3665, -1.3323, -1.3068, -1.2274, -1.2425, -1.2308, -1.3553, -1.3102, -1.232, -1.4006, -1.6049, -1.6494, -1.7022, -1.6922, -1.759, -1.7027, -2.2783, -3.9462, -5.8421, -8.1127, -10.6619, -13.6482, -16.8177, -20.6102, -25.8113, -32.0407, -40.2741, -51.0496, -66.2185, -90.4265, -125.35, -189.042};
static const double AU_IM[49] = {-3.0413, -3.1759, -3.2857, -3.3277, -3.3968, -3.3904, -3.51, -3.6062, -3.7102, -3.8252, -3.8922, -4.055, -4.1633, -4.3385, -4.4901, -4.7223, -4.9763, -5.2824, -5.4949, -5.5964, -5.7803, -5.7926, -5.8458, -5.5737, -5.5382, -5.598, -5.6092, -5.6444, -5.7389, -5.7174, -5.6492, -5.2826, -4.8444, -3.8126, -2.5804, -2.1113, -1.6605, -1.3742, -1.0352, -1.0668, -1.2718, -1.6266, -1.9254, -2.794, -3.861, -5.7015, -8.1863, -12.5552, -25.3552};

static const double CU_XI[49] = {187.9, 191.6, 195.3, 199.3, 203.3, 207.3, 211.9, 216.4, 221.4, 226.2, 231.3, 237.1, 242.6, 249, 255.1, 261.6, 268.9, 276.1, 284.4, 292.4, 300.9, 310.7, 320.4, 331.5, 342.5, 354.2, 367.9, 381.5, 397.4, 413.3, 430.5, 450.9, 471.4, 495.9, 520.9, 548.6, 582.1, 616.8, 659.5, 704.5, 756, 821.1, 892, 984, 1088, 1216, 1393, 1610, 1937};
static const double CU_RE[49] = {-0.904, -1.024, -1.1327, -1.2686, -1.4224, -1.5367, -1.6442, -1.7202, -1.7403, -1.7334, -1.6984, -1.6088, -1.4408, -1.3022, -1.043, -0.8714, -0.6797, -0.5777, -0.5642, -0.6503, -0.859, -1.085, -1.2747, -1.5204, -1.6249, -1.7942, -2.051, -2.4131, -2.7351, -3.2324, -3.7505, -4.208, -4.6028, -5.0857, -5.4093, -5.6005, -6.8216, -10.182, -13.9916, -17.6379, -21.7046, -26.7648, -33.1798, -41.1268, -51.9555, -67.7496, -88.7347, -123.077, -179.177};
static const double CU_IM[49] = {-2.5136, -2.6372, -2.7936, -2.9263, -3.069, -3.23, -3.4341, -3.6698, -3.9256, -4.1725, -4.4083, -4.6131, -4.8213, -4.8854, -4.9096, -4.7686, -4.8372, -4.8063, -4.7357, -4.6377, -4.7012, -4.772, -4.9211, -4.8803, -5.0701, -5.2498, -5.372, -5.4397, -5.5862, -5.6499, -5.7625, -5.9446, -6.2075, -6.2562, -6.1549, -5.2571, -3.7856, -1.923, -1.6487, -1.7661, -2.2392, -2.6936, -3.4608, -4.1094, -5.1962, -7.9152, -11.3268, -16.9024, -29.2774};

/* локальная линейная интерполяция Лагранжа (mint=2), порт из EPS_MET_AuAgCu.c */
static double lagr(double x, const double *xi, const double *col, int n)
{
	int i = 0, b1, ed, ii, jj, q = n - 1;
	double total = 0.0, term;
	while (i < n && xi[i] < x) i++;
	if (i >= n) i = n - 1;
	b1 = i - 1;               /* mint/2 = 1 */
	if (b1 < 0) b1 = 0;
	if (b1 + 1 > q) b1 = q - 1;
	ed = b1 + 1;
	for (ii = b1; ii <= ed; ii++) {
		term = col[ii];
		for (jj = b1; jj <= ed; jj++)
			if (jj != ii) term *= (x - xi[jj]) / (xi[ii] - xi[jj]);
		total += term;
	}
	return total;
}

/* диспетчер материала: n_mat = 4 (Au), 1 (Cu) — минимальная демо-версия */
void libmat(int n_mat, double lam, double *Re_eps, double *Im_eps)
{
	if (n_mat == 4) {          /* Au */
		*Re_eps = lagr(lam, AU_XI, AU_RE, 49);
		*Im_eps = lagr(lam, AU_XI, AU_IM, 49);
	} else if (n_mat == 1) {   /* Cu */
		*Re_eps = lagr(lam, CU_XI, CU_RE, 49);
		*Im_eps = lagr(lam, CU_XI, CU_IM, 49);
	} else {
		printf("libmat(demo): n_mat=%d не поддержан (только 4=Au, 1=Cu)\n", n_mat);
		*Re_eps = *Im_eps = 0.0;
	}
}

/* --- собираем реальные функции моделей --- */
#include "../source_c/COMPOSIT.c"
#include "../source_c/COMPOSITE_MODELS.c"

/* поиск пика |Im eps_eff| одной моделью на сетке длин волн */
typedef void (*model1)(double, double, double, int, double *, double *);

static void sweep1(const char *label, model1 fn, double emR, double emI,
	double C, int n_mat, double lo, double hi, double step)
{
	double peak_lam = lo, peak_absIm = -1.0, peak_Im = 0.0, lam, ocR, ocI;
	for (lam = lo; lam <= hi + 1e-9; lam += step) {
		f = lam;                       /* обёртки читают глобальную f */
		fn(emR, emI, C, n_mat, &ocR, &ocI);
		if (fabs(ocI) > peak_absIm) { peak_absIm = fabs(ocI); peak_Im = ocI; peak_lam = lam; }
	}
	printf("  %-26s пик %.0f нм,  Im eps_eff = %+.3f\n", label, peak_lam, peak_Im);
}

/* MLWA имеет доп. параметр R_p — своя обёртка sweep */
static void sweep_mlwa(double R_p, double emR, double emI, double C, int n_mat,
	double lo, double hi, double step)
{
	double peak_lam = lo, peak_absIm = -1.0, peak_Im = 0.0, lam, ocR, ocI;
	char label[48];
	for (lam = lo; lam <= hi + 1e-9; lam += step) {
		f = lam;
		COMPOSITE_MLWA(emR, emI, C, n_mat, R_p, &ocR, &ocI);
		if (fabs(ocI) > peak_absIm) { peak_absIm = fabs(ocI); peak_Im = ocI; peak_lam = lam; }
	}
	sprintf(label, "MLWA, R_p=%.0f нм", R_p);
	printf("  %-26s пик %.0f нм,  Im eps_eff = %+.3f\n", label, peak_lam, peak_Im);
}

int main(void)
{
	double n_host = 1.77;
	double emR = n_host * n_host, emI = 0.0;   /* eps матрицы = 3.1329 */
	double C = 0.10;                           /* объёмная доля частиц */
	double ocR, ocI;
	int mats[2] = {4, 1};
	const char *names[2] = {"Au", "Cu"};
	int mi;

	printf("Пример: eps_eff композита (n_host=%.2f, C=%.0f%%), 3 модели эффективной среды\n", n_host, C * 100);
	printf("Материалы: Johnson & Christy 1972 (Au, Cu). Сравнение с таблицей 2 статьи.\n\n");

	/* --- одиночный вызов каждой модели (пример использования) --- */
	f = 535.0;
	printf("Точечный пример на 535 нм, Au (n_mat=4):\n");
	COMPOSITE(emR, emI, C, 4, &ocR, &ocI);
	printf("  COMPOSITE      -> eps_eff = %+.4f %+.4fi\n", ocR, ocI);
	BRUGGEMAN(emR, emI, C, 4, &ocR, &ocI);
	printf("  BRUGGEMAN      -> eps_eff = %+.4f %+.4fi\n", ocR, ocI);
	COMPOSITE_MLWA(emR, emI, C, 4, 30.0, &ocR, &ocI);
	printf("  COMPOSITE_MLWA -> eps_eff = %+.4f %+.4fi   (R_p=30 нм)\n\n", ocR, ocI);

	/* --- воспроизведение таблицы 2: пик |Im eps_eff| по каждой модели --- */
	for (mi = 0; mi < 2; mi++) {
		printf("%s (n_mat=%d):\n", names[mi], mats[mi]);
		sweep1("Клаузиус-Моссотти (1)", COMPOSITE,  emR, emI, C, mats[mi], 300.0, 900.0, 2.0);
		sweep1("Бруггеман (4)",         BRUGGEMAN,  emR, emI, C, mats[mi], 300.0, 900.0, 2.0);
		sweep_mlwa(10.0, emR, emI, C, mats[mi], 300.0, 900.0, 2.0);
		sweep_mlwa(30.0, emR, emI, C, mats[mi], 300.0, 900.0, 2.0);
		sweep_mlwa(50.0, emR, emI, C, mats[mi], 300.0, 900.0, 2.0);
		printf("\n");
	}

	printf("Ожидание (таблица 2 статьи):\n");
	printf("  Au: КМ 572/-6.05, Бруггеман 708/-2.66, MLWA 580/-6.69, 640/-6.49, 790/-3.59\n");
	printf("  Cu: КМ 596/-3.32, Бруггеман 702/-2.56, MLWA 600/-3.80, 640/-5.60, 786/-3.37\n");
	return 0;
}
