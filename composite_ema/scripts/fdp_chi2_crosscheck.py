"""Сверка записи Лерера через chi_2 с формулой Максвелла-Гарнетта его же кода.

Запись из «Композит из наночастиц. Диэлектрическая проницаемость» (А.М. Лерер):

    eps_c = 1 + (1 - C)(eps_m - 1) + C * chi_2,
    chi_2 = 3 eps_m (eps_n - eps_m)/(eps_n + 2 eps_m).

Код COMPOSIT.c реализует Максвелла-Гарнетта:

    eps_eff = eps_m (1 + 2 C K)/(1 - C K),  K = (eps_n - eps_m)/(eps_n + 2 eps_m).

Скрипт проверяет три вещи:
 1) тривиальный контроль eps_n = eps_m (частицы неотличимы от матрицы);
 2) предел малых C у обеих записей;
 3) численное расхождение при рабочих параметрах статьи.
"""

from __future__ import annotations

import composite_ema as ce


def mg(eps_m: complex, eps_n: complex, C: float) -> complex:
    K = (eps_n - eps_m) / (eps_n + 2.0 * eps_m)
    return eps_m * (1.0 + 2.0 * C * K) / (1.0 - C * K)


def chi2(eps_m: complex, eps_n: complex) -> complex:
    return 3.0 * eps_m * (eps_n - eps_m) / (eps_n + 2.0 * eps_m)


def lerer_chi2_form(eps_m: complex, eps_n: complex, C: float) -> complex:
    """Запись из записки Лерера, буквально."""
    return 1.0 + (1.0 - C) * (eps_m - 1.0) + C * chi2(eps_m, eps_n)


def dilute_form(eps_m: complex, eps_n: complex, C: float) -> complex:
    """Разложение Максвелла-Гарнетта до первого порядка по C."""
    return eps_m + C * chi2(eps_m, eps_n)


def main() -> None:
    eps_m = complex(1.77**2, 0.0)

    print("=" * 72)
    print("1) Контроль: частицы неотличимы от матрицы, eps_n = eps_m")
    print("   Физически обязан получиться ровно eps_m = %.4f" % eps_m.real)
    print("=" * 72)
    for C in (0.01, 0.10, 0.30):
        a = mg(eps_m, eps_m, C)
        b = lerer_chi2_form(eps_m, eps_m, C)
        c = dilute_form(eps_m, eps_m, C)
        print(
            "   C=%.2f:  Максвелл-Гарнетт=%.6f   запись chi_2=%.6f   eps_m + C*chi_2=%.6f"
            % (C, a.real, b.real, c.real)
        )
    print("   -> у записи chi_2 появляется ложный сдвиг -C(eps_m - 1) = %.4f при C=0,1"
          % (-0.10 * (eps_m.real - 1.0)))

    print()
    print("=" * 72)
    print("2) Предел C -> 0: производная d eps / dC")
    print("=" * 72)
    eps_n = ce.eps_au(535.0)
    h = 1e-7
    d_mg = (mg(eps_m, eps_n, h) - eps_m) / h
    d_lerer = (lerer_chi2_form(eps_m, eps_n, h) - lerer_chi2_form(eps_m, eps_n, 0.0)) / h
    print("   Максвелл-Гарнетт : %+.6f%+.6fj" % (d_mg.real, d_mg.imag))
    print("   chi_2 (Лерер)    : %+.6f%+.6fj" % (d_lerer.real, d_lerer.imag))
    print("   chi_2 сам по себе: %+.6f%+.6fj" % (chi2(eps_m, eps_n).real, chi2(eps_m, eps_n).imag))
    print("   разность         : %+.6f%+.6fj  = -(eps_m - 1) = %+.6f"
          % ((d_lerer - d_mg).real, (d_lerer - d_mg).imag, -(eps_m.real - 1.0)))

    print()
    print("=" * 72)
    print("3) Рабочая точка статьи: n_h=1,77, C=0,10, Au, 300-900 нм")
    print("=" * 72)
    print("   lam, нм |      Максвелл-Гарнетт |         запись chi_2 |  |разность|")
    worst = 0.0
    for lam in (400.0, 500.0, 535.0, 572.0, 600.0, 700.0, 800.0, 900.0):
        eps_n = ce.eps_au(lam)
        a = mg(eps_m, eps_n, 0.10)
        b = lerer_chi2_form(eps_m, eps_n, 0.10)
        worst = max(worst, abs(a - b))
        print(
            "   %7.0f | %+8.4f%+8.4fj | %+8.4f%+8.4fj | %8.4f"
            % (lam, a.real, a.imag, b.real, b.imag, abs(a - b))
        )
    print("   максимальное расхождение по этим точкам: %.4f" % worst)

    print()
    print("=" * 72)
    print("4) Правильная малая-C запись: eps_m + C*chi_2 (без множителя (1-C) у матрицы)")
    print("=" * 72)
    for C in (0.01, 0.05, 0.10):
        eps_n = ce.eps_au(535.0)
        a = mg(eps_m, eps_n, C)
        c = dilute_form(eps_m, eps_n, C)
        print(
            "   C=%.2f:  Максвелл-Гарнетт=%+.5f%+.5fj   eps_m+C*chi_2=%+.5f%+.5fj   |разность|=%.5f"
            % (C, a.real, a.imag, c.real, c.imag, abs(a - c))
        )
    print("   -> совпадение по первому порядку C, расхождение растёт как C^2 (это ожидаемо).")


if __name__ == "__main__":
    main()
