/*
  Дополнительные модели эффективной среды для слоя диэлектрика с плазмонными
  наночастицами (ПНЧ). Расширение COMPOSIT.c: рядом с формулой Клаузиуса-Моссотти
  (Максвелла-Гарнетта) — COMPOSITE()/COMPOSITE_3() — добавлены две модели,
  использованные в статье для оценки модельной неопределённости:

    BRUGGEMAN       — симметричная модель Бруггемана, формула (4) статьи;
    COMPOSITE_MLWA  — Максвелл-Гарнетт с размерной поправкой в модифицированном
                      длинноволновом приближении (MLWA): динамическая
                      деполяризация + радиационное затухание, ссылка [9] статьи.

  Статья: А.М. Лерер, А.Б. Клещенков, Н.В. Горбачев, «Влияние модели эффективной
  среды и источника оптических констант на расчётные характеристики поглощающих
  метаповерхностей с плазмонными наночастицами» (Оптика и спектроскопия).

  Соглашения — те же, что в COMPOSIT.c:
   - временная зависимость e^{+i w t}; у пассивной среды eps = eps' - i*eps'',
     Im(eps) <= 0 (мнимая часть, которую возвращает libmat(), отрицательна);
   - divide_(x,y,x1,y1,&u,&v) = (x + i y)/(x1 + i y1);
   - mult_  (x,y,x1,y1,&u,&v) = (x + i y)*(x1 + i y1);
   - libmat(n_mat, f, &Re_eps,&Im_eps) — выбор материала частиц из базы;
     глобальная f — длина волны в нм (как в вызове из COMPOSITE()).

  Для каждой модели есть две функции:
   - *_eps — «ядро», принимает eps частиц напрямую (emR,emI — матрица;
     epR,epI — частицы); удобно для автономной проверки без базы;
   - обёртка BRUGGEMAN()/COMPOSITE_MLWA() — берёт eps частиц из базы через
     libmat(n_mat, f, ...), сигнатура согласована с COMPOSITE().

  Численно проверено против независимого Python-порта
  (composite_ema/scripts/compare_effective_models.py): совпадение до 6 знаков
  на контрольных точках, воспроизводит таблицу 2 статьи (положения и амплитуды
  пиков |Im eps_eff| для Au и Cu; Бруггеман 708/702 нм, MLWA R=10/30/50 нм —
  580/640/790 нм и т.д.).
*/

#include "VIBRr.h"
#include <math.h>

/* Главное значение комплексного корня sqrt(a + i b). */
static void csqrt_(double a, double b, double *sr, double *si)
{
	double r = sqrt(a * a + b * b);
	double g = sqrt(0.5 * (r + a));
	double d = sqrt(0.5 * (r - a));
	if (b < 0.0) d = -d;
	*sr = g;
	*si = d;
}

/*
  Симметричная модель Бруггемана (формула (4) статьи):
     C*(eps_p - eps)/(eps_p + 2 eps) + (1-C)*(eps_m - eps)/(eps_m + 2 eps) = 0.
  Приводится к квадратному уравнению
     2 eps^2 - b eps - eps_p*eps_m = 0,  b = (3C-1) eps_p + (2-3C) eps_m,
  корни  eps = ( b +- sqrt(b^2 + 8 eps_p eps_m) ) / 4.
  Физическая ветвь — корень, ближайший к проницаемости матрицы eps_m
  (для малых C непрерывно переходит в eps_m).
*/
void BRUGGEMAN_eps(double emR, double emI,   // eps матрицы
	double epR, double epI,                  // eps материала частиц
	double C,                                // объёмная доля частиц
	double *ocR, double *ocI)                // eps композита
{
	double a1 = 3.0 * C - 1.0, a2 = 2.0 - 3.0 * C;
	double bR = a1 * epR + a2 * emR;
	double bI = a1 * epI + a2 * emI;
	double b2R, b2I, pmR, pmI, discR, discI, sR, sI;
	double r1R, r1I, r2R, r2I, d1, d2, dr, di;

	mult_(bR, bI, bR, bI, &b2R, &b2I);        // b^2
	mult_(epR, epI, emR, emI, &pmR, &pmI);    // eps_p * eps_m
	discR = b2R + 8.0 * pmR;
	discI = b2I + 8.0 * pmI;
	csqrt_(discR, discI, &sR, &sI);

	r1R = 0.25 * (bR + sR); r1I = 0.25 * (bI + sI);
	r2R = 0.25 * (bR - sR); r2I = 0.25 * (bI - sI);

	dr = r1R - emR; di = r1I - emI; d1 = dr * dr + di * di;
	dr = r2R - emR; di = r2I - emI; d2 = dr * dr + di * di;

	if (d1 <= d2) { *ocR = r1R; *ocI = r1I; }
	else          { *ocR = r2R; *ocI = r2I; }
}//BRUGGEMAN_eps

void BRUGGEMAN(double R_eps_m, double I_eps_m, //Комплексная диэл. проницаемость матрицы
	double C,                                  //концентрация н/частиц
	int n_mat,                                 //номер материала н/частиц в базе (см. libmat)
	double *R_eps_c, double *I_eps_c)          //Комплексная диэл. проницаемость композита
{
	double Re_eps, Im_eps;
	libmat(n_mat, f, &Re_eps, &Im_eps);
	BRUGGEMAN_eps(R_eps_m, I_eps_m, Re_eps, Im_eps, C, R_eps_c, I_eps_c);
}//BRUGGEMAN


/*
  Максвелл-Гарнетт с размерной поправкой (MLWA) для сферической частицы радиуса
  R_p. Квазистатическая поляризуемость заменяется на MLWA-поляризуемость:
     eta = C*K / [ 1 - x^2 K - i (2/3) x^3 K ],  K = (eps_p - eps_m)/(eps_p + 2 eps_m),
  где x = k_h R_p = 2 pi n_h R_p / lam — безразмерный размерный параметр в матрице,
  n_h = Re( sqrt(eps_m) ). Эффективная проницаемость: eps_eff = eps_m (1 + 2 eta)/(1 - eta).
  Формула MLWA записана в конвенции e^{-i w t} (у пассивной среды Im>0), поэтому
  внутри функции знак Im инвертируется на входе и возвращается обратно на выходе.
  R_p и lam — в одинаковых единицах (нм); в отношении R_p/lam единицы сокращаются.
*/
void COMPOSITE_MLWA_eps(double emR, double emI,  // eps матрицы (конвенция e^{+iwt})
	double epR, double epI,                      // eps материала частиц (конвенция e^{+iwt})
	double C,                                     // объёмная доля частиц
	double R_p,                                   // радиус частицы, нм
	double lam,                                   // длина волны, нм
	double *ocR, double *ocI)                     // eps композита (конвенция e^{+iwt})
{
	double hR = emR, hI = -emI;   // eps_m в конвенции e^{-i w t}
	double pR = epR, pI = -epI;   // eps_p в конвенции e^{-i w t}
	double sr, si, n_h, x, x2, x3;
	double numR, numI, denR, denI, KR, KI, dR, dI, etaR, etaI, rR, rI, ecR, ecI;
	double pi = 3.14159265358979323846;

	csqrt_(hR, hI, &sr, &si);     // n_h = Re( sqrt(eps_m) )
	n_h = sr;
	x = 2.0 * pi * n_h * R_p / lam;
	x2 = x * x; x3 = x2 * x;

	// K = (eps_p - eps_m)/(eps_p + 2 eps_m)
	numR = pR - hR;        numI = pI - hI;
	denR = pR + 2.0 * hR;  denI = pI + 2.0 * hI;
	divide_(numR, numI, denR, denI, &KR, &KI);

	// знаменатель MLWA: 1 - x^2 K - i (2/3) x^3 K
	dR = 1.0 - x2 * KR + (2.0 / 3.0) * x3 * KI;
	dI =     - x2 * KI - (2.0 / 3.0) * x3 * KR;

	// eta = C*K / denominator
	divide_(C * KR, C * KI, dR, dI, &etaR, &etaI);

	// eps_eff = eps_m (1 + 2 eta)/(1 - eta)
	numR = 1.0 + 2.0 * etaR; numI = 2.0 * etaI;
	denR = 1.0 - etaR;       denI = -etaI;
	divide_(numR, numI, denR, denI, &rR, &rI);
	mult_(hR, hI, rR, rI, &ecR, &ecI);

	*ocR = ecR;
	*ocI = -ecI;                  // обратно в конвенцию e^{+i w t}
}//COMPOSITE_MLWA_eps

void COMPOSITE_MLWA(double R_eps_m, double I_eps_m, //Комплексная диэл. проницаемость матрицы
	double C,                                       //концентрация н/частиц
	int n_mat,                                      //номер материала н/частиц в базе (см. libmat)
	double R_p,                                     //радиус частицы, нм
	double *R_eps_c, double *I_eps_c)               //Комплексная диэл. проницаемость композита
{
	double Re_eps, Im_eps;
	libmat(n_mat, f, &Re_eps, &Im_eps);
	COMPOSITE_MLWA_eps(R_eps_m, I_eps_m, Re_eps, Im_eps, C, R_p, f, R_eps_c, I_eps_c);
}//COMPOSITE_MLWA
