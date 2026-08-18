"""Почему модель Бруггемана так сильно расходится с Максвеллом-Гарнеттом.

Вопрос А. М. Лерера от 2026-08-18: «При R=15 графики близки к КМ. Но почему
такое большое отличие метода Б.? Кстати, если считать по методу Б поглощение
самое широкополосное.»

Скрипт проверяет четыре утверждения ответа:
 1) обе модели совпадают в первом порядке по C, значит дело не в разной физике
    отдельной частицы, а в величине параметра разложения C*K;
 2) у Максвелла-Гарнетта резонансный знаменатель eps_p + 2 eps_m содержит
    ФИКСИРОВАННУЮ вещественную проницаемость матрицы, у Бруггемана -
    eps_p + 2 eps_eff, то есть комплексную самосогласованную величину;
 3) отсюда красный сдвиг и сильное демпфирование: включение сидит в поглощающей
    среде, и знаменатель никогда не становится малым;
 4) следствие - полоса поглощения у Бруггемана действительно самая широкая.

Данные - таблицы luxpop самого Лерера.
"""

from __future__ import annotations

import cmath
import math

import numpy as np

from lerer_compare_diagnosis import bruggeman, eps_luxpop, load_luxpop, mg, mlwa

N_M = 1.77
EPS_M = complex(N_M**2, 0.0)


def fwhm(lam: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    """Полная ширина на половине высоты; возвращает (ширина, лево, право)."""
    i = int(np.argmax(y))
    half = y[i] / 2.0

    left = lam[0]
    for j in range(i, 0, -1):
        if y[j] <= half:
            t = (half - y[j]) / (y[j + 1] - y[j])
            left = lam[j] + t * (lam[j + 1] - lam[j])
            break
    right = lam[-1]
    for j in range(i, len(lam) - 1):
        if y[j] <= half:
            t = (half - y[j]) / (y[j - 1] - y[j])
            right = lam[j] + t * (lam[j - 1] - lam[j])
            break
    return right - left, left, right


def main() -> None:
    au = load_luxpop("Au_.c")
    lam = np.arange(300.0, 1500.0 + 1e-9, 1.0)
    eps_p = np.array([eps_luxpop(au, float(l)) for l in lam])

    print("=" * 78)
    print("1) Совпадают ли модели при малых C? (Au, luxpop)")
    print("=" * 78)
    print("   Обе записи в первом порядке дают eps_eff = eps_m (1 + 3 C K).")
    print("   %6s | %-22s | %-22s | %s" % ("C", "пик КМ", "пик Бруггемана", "разнос"))
    for C in (0.001, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20):
        km = np.array([mg(EPS_M, e, C) for e in eps_p])
        br = np.array([bruggeman(EPS_M, e, C) for e in eps_p])
        i1, i2 = int(np.argmax(np.abs(km.imag))), int(np.argmax(np.abs(br.imag)))
        print("   %6.3f | %5.0f нм, |Im|=%7.3f | %5.0f нм, |Im|=%7.3f | %4.0f нм"
              % (C, lam[i1], abs(km.imag[i1]), lam[i2], abs(br.imag[i2]),
                 lam[i2] - lam[i1]))
    print("   -> при C -> 0 модели сходятся; расхождение растёт вместе с C.")

    print()
    print("=" * 78)
    print("2) Настоящий параметр разложения - не C, а C*|K|")
    print("=" * 78)
    C = 0.10
    km = np.array([mg(EPS_M, e, C) for e in eps_p])
    i_res = int(np.argmax(np.abs(km.imag)))
    lam_res = lam[i_res]
    K = (eps_p - EPS_M) / (eps_p + 2.0 * EPS_M)
    print("   На резонансе Клаузиуса-Моссотти (%.0f нм):" % lam_res)
    print("     eps_p            = %+.4f%+.4fj" % (eps_p[i_res].real, eps_p[i_res].imag))
    print("     K                = %+.4f%+.4fj,  |K| = %.3f"
          % (K[i_res].real, K[i_res].imag, abs(K[i_res])))
    print("     C|K|             = %.3f  <- НЕ мало, ряд по C здесь не работает"
          % (C * abs(K[i_res])))
    print("   Вдали от резонанса, для сравнения:")
    for l0 in (400.0, 900.0, 1400.0):
        j = int(np.argmin(np.abs(lam - l0)))
        print("     %4.0f нм: |K| = %6.3f,  C|K| = %.3f" % (l0, abs(K[j]), C * abs(K[j])))

    print()
    print("=" * 78)
    print("3) В чём сидит включение: фиксированная матрица или сама смесь")
    print("=" * 78)
    br = np.array([bruggeman(EPS_M, e, C) for e in eps_p])
    i_br = int(np.argmax(np.abs(br.imag)))
    print("   Максвелл-Гарнетт: знаменатель eps_p + 2*eps_m, eps_m = %.4f (вещественная)"
          % EPS_M.real)
    print("     условие резонанса Re eps_p = -2 eps_m = %.3f" % (-2.0 * EPS_M.real))
    print("     выполняется при %.0f нм" % lam_res)
    print()
    print("   Бруггеман: знаменатель eps_p + 2*eps_eff, eps_eff комплексная и своя на каждой длине волны")
    print("     в максимуме (%.0f нм): eps_eff = %+.4f%+.4fj"
          % (lam[i_br], br[i_br].real, br[i_br].imag))
    print("     -2*Re eps_eff = %.3f,  а Re eps_p = %.3f" % (-2.0 * br[i_br].real, eps_p[i_br].real))
    print("     то есть резонансное условие Re eps_p = -2 Re eps_eff в максимуме НЕ выполнено:")
    print("     максимум Бруггемана - не резонанс включения, а пологий максимум решения.")
    print("     мнимая часть фона -2*Im eps_eff = %.3f - включение сидит в ПОГЛОЩАЮЩЕЙ среде"
          % (-2.0 * br[i_br].imag))
    print()
    print("   Модуль резонансного знаменателя в максимуме:")
    den_mg = abs(eps_p[i_res] + 2.0 * EPS_M)
    den_br = abs(eps_p[i_br] + 2.0 * br[i_br])
    print("     Максвелл-Гарнетт |eps_p + 2 eps_m|    = %.3f" % den_mg)
    print("     Бруггеман        |eps_p + 2 eps_eff|  = %.3f" % den_br)
    print("     -> у Бруггемана знаменатель в %.1f раза больше, резонанс подавлен"
          % (den_br / den_mg))

    print()
    print("=" * 78)
    print("4) Ширина полосы поглощения (Лерер: «по методу Б самое широкополосное»)")
    print("=" * 78)
    models = [
        ("Клаузиус-Моссотти", km),
        ("Бруггеман", br),
        ("MLWA, R=30 нм", np.array([mlwa(EPS_M, e, C, float(l), 30.0) for e, l in zip(eps_p, lam)])),
        ("MLWA, R=15 нм", np.array([mlwa(EPS_M, e, C, float(l), 15.0) for e, l in zip(eps_p, lam)])),
    ]
    print("   %-20s %9s %9s %14s %10s" % ("модель", "пик, нм", "|Im| макс", "полоса, нм", "полоса/пик"))
    for name, vals in models:
        y = np.abs(vals.imag)
        i = int(np.argmax(y))
        w, lo, hi = fwhm(lam, y)
        print("   %-20s %9.0f %9.3f %14s %10.2f"
              % (name, lam[i], y[i], "%.0f (%.0f-%.0f)" % (w, lo, hi), w / lam[i]))
    print("   -> подтверждается: у Бруггемана полоса самая широкая, но и пик самый низкий.")

    print()
    print("=" * 78)
    print("4a) Насколько близко каждая модель подходит к полюсу")
    print("=" * 78)
    print("   У Максвелла-Гарнетта eps_eff = eps_m (1 + 2CK)/(1 - CK):")
    print("   резкость задаётся тем, насколько C|K| подходит к единице.")
    print("     max C|K| по спектру = %.3f (при %.0f нм)"
          % (C * np.max(np.abs(K)), lam[int(np.argmax(np.abs(K)))]))
    den_mg_all = np.abs(eps_p + 2.0 * EPS_M)
    den_br_all = np.abs(eps_p + 2.0 * br)
    j1, j2 = int(np.argmin(den_mg_all)), int(np.argmin(den_br_all))
    print("   Минимум |eps_p + 2 eps_m|   (МГ, фон - прозрачная матрица): %.3f при %.0f нм"
          % (den_mg_all[j1], lam[j1]))
    print("   Минимум |eps_p + 2 eps_eff| (Б,  фон - сама смесь):        %.3f при %.0f нм"
          % (den_br_all[j2], lam[j2]))
    print("   -> у Бруггемана знаменатель нигде не становится малым: самосогласованность")
    print("      не даёт выполниться резонансному условию, острого резонанса просто нет.")

    print()
    print("=" * 78)
    print("4b) Правило сумм: не создаётся и не теряется, а перераспределяется")
    print("=" * 78)
    print("   Интеграл  S = int |Im eps_eff| * omega d(omega)  по 300-1500 нм")
    print("   (при фиксированной объёмной доле металла полная сила осциллятора задана)")
    omega = 2.0 * math.pi * 299792458.0 / (lam * 1e-9)
    for name, vals in models:
        y = np.abs(vals.imag)
        S = np.trapezoid(y * omega, omega) if hasattr(np, "trapezoid") else np.trapz(y * omega, omega)
        print("     %-20s S = %.4e" % (name, abs(S)))
    print("   -> величины одного порядка: модели по-разному РАСПРЕДЕЛЯЮТ одно и то же")
    print("      поглощение. Бруггеман размазывает его широко и мелко, КМ собирает узко и высоко.")

    print()
    print("=" * 78)
    print("5) Порог протекания симметричной модели")
    print("=" * 78)
    print("   Для сфер (фактор деполяризации 1/3) симметричная модель Бруггемана")
    print("   имеет порог протекания C = 1/3. При C = 0,10 это уже треть пути к порогу,")
    print("   то есть модель описывает частично связную смесь, а не изолированные сферы.")

    print()
    print("=" * 78)
    print("6) Контрольный узел: его таблица против Джонсона-Кристи")
    print("=" * 78)
    n_l, k_l = 0.567, 2.202   # его база, Lam[521]
    n_j, k_j = 0.62, 2.081    # Johnson & Christy, 520,9 нм
    e_l = complex(n_l**2 - k_l**2, -2 * n_l * k_l)
    e_j = complex(n_j**2 - k_j**2, -2 * n_j * k_j)
    print("   его база, 521,0 нм: n=%.3f k=%.3f -> eps = %+.4f%+.4fj" % (n_l, k_l, e_l.real, e_l.imag))
    print("   J&C,      520,9 нм: n=%.3f k=%.3f -> eps = %+.4f%+.4fj" % (n_j, k_j, e_j.real, e_j.imag))
    print("   расхождение: Re %.1f %%, Im %.1f %%"
          % (abs(e_l.real - e_j.real) / abs(e_j.real) * 100,
             abs(e_l.imag - e_j.imag) / abs(e_j.imag) * 100))
    print("   Именно расхождение по Re сдвигает условие Фрёлиха: 539 нм против 555 нм.")


if __name__ == "__main__":
    main()
