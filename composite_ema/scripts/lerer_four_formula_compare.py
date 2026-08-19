"""Сравнение четырёх формул эффективной среды по письму А.М. Лерера 2026-08-19.

Повод: файлы «Сравнение 4 формул расчета композита.docx» и
«Композит из наночастиц_Диэлектрическая проницаемость.docx». Лерер сопоставил
на своей метаповерхности четыре кривые поглощения P(lambda): «Клаузиус-Моссотти
(КМ)», «Бруггемана (Б)», «Максвелл-Гарнетт» и «Лерер».

Задача скрипта - не воспроизвести RCWA-расчёт метаповерхности, а установить
точные соотношения между самими формулами eps_eff при параметрах статьи
(n_h = 1.77, C = 0.10, Au по Johnson-Christy), чтобы можно было судить, какие
пары кривых обязаны совпадать, а какие - расходиться.

Конвенция вывода: exp(-i w t), у пассивной среды Im eps >= 0. Данные
composite_ema хранятся в конвенции Лерера exp(+i w t), поэтому сопрягаются.

Формулы:
  MG      - полный Максвелл-Гарнетт: eps_h (1 + 2CK)/(1 - CK)
  DILUTE  - первый порядок по C:     eps_h + C * chi2,  chi2 = 3 eps_h K
  LERER   - запись Лерера:           1 + (1-C)(eps_h - 1) + C * chi2
  BRUG    - симметричный Бруггеман
  LL      - Лоренц-Лоренц / классический Клаузиус-Моссотти относительно вакуума:
            (eps-1)/(eps+2) = (1-C)(eps_h-1)/(eps_h+2) + C (eps_p-1)/(eps_p+2)
"""

from __future__ import annotations

import cmath

import composite_ema as ce

N_H = 1.77
EPS_H = complex(N_H**2, 0.0)
C = 0.10
LAM_LO, LAM_HI, LAM_STEP = 400.0, 750.0, 1.0


def eps_au(lam_nm: float) -> complex:
    """Au по Johnson-Christy в конвенции exp(-i w t): Im eps >= 0."""
    return ce.eps_au(lam_nm).conjugate()


def k_contrast(eps_h: complex, eps_p: complex) -> complex:
    return (eps_p - eps_h) / (eps_p + 2.0 * eps_h)


def chi2(eps_h: complex, eps_p: complex) -> complex:
    return 3.0 * eps_h * k_contrast(eps_h, eps_p)


def mg(eps_h: complex, eps_p: complex, c: float) -> complex:
    K = k_contrast(eps_h, eps_p)
    return eps_h * (1.0 + 2.0 * c * K) / (1.0 - c * K)


def dilute(eps_h: complex, eps_p: complex, c: float) -> complex:
    return eps_h + c * chi2(eps_h, eps_p)


def lerer(eps_h: complex, eps_p: complex, c: float) -> complex:
    return 1.0 + (1.0 - c) * (eps_h - 1.0) + c * chi2(eps_h, eps_p)


def bruggeman(eps_h: complex, eps_p: complex, c: float) -> complex:
    b = (3.0 * c - 1.0) * eps_p + (2.0 - 3.0 * c) * eps_h
    disc = cmath.sqrt(b * b + 8.0 * eps_p * eps_h)
    r1, r2 = (b + disc) / 4.0, (b - disc) / 4.0
    return r1 if abs(r1 - eps_h) <= abs(r2 - eps_h) else r2


def lorentz_lorenz(eps_h: complex, eps_p: complex, c: float) -> complex:
    """Классический Клаузиус-Моссотти относительно вакуума."""
    s = (1.0 - c) * (eps_h - 1.0) / (eps_h + 2.0) + c * (eps_p - 1.0) / (eps_p + 2.0)
    return (1.0 + 2.0 * s) / (1.0 - s)


MODELS = {
    "MG": mg,
    "DILUTE": dilute,
    "LERER": lerer,
    "BRUG": bruggeman,
    "LL": lorentz_lorenz,
}


def grid() -> list[float]:
    n = int(round((LAM_HI - LAM_LO) / LAM_STEP)) + 1
    return [LAM_LO + i * LAM_STEP for i in range(n)]


def main() -> None:
    lams = grid()
    eps_p = {lam: eps_au(lam) for lam in lams}
    vals = {name: [fn(EPS_H, eps_p[lam], C) for lam in lams] for name, fn in MODELS.items()}

    print("=" * 78)
    print("Параметры: n_h = %.2f, eps_h = %.4f, C = %.2f, Au (Johnson-Christy)" % (N_H, EPS_H.real, C))
    print("Диапазон %.0f-%.0f нм, шаг %.0f нм. Конвенция exp(-i w t), Im eps >= 0." % (LAM_LO, LAM_HI, LAM_STEP))
    print("=" * 78)

    # 1. Тождество LERER = DILUTE - C (eps_h - 1)
    offset = C * (EPS_H.real - 1.0)
    worst = max(abs(l - (d - offset)) for l, d in zip(vals["LERER"], vals["DILUTE"]))
    print("\n1) Тождество  LERER = DILUTE - C(eps_h - 1)")
    print("   сдвиг C(eps_h - 1) = %.6f (действительный, не зависит от длины волны)" % offset)
    print("   max |LERER - (DILUTE - сдвиг)| по всей сетке = %.3e" % worst)
    print("   => Im LERER == Im DILUTE тождественно; расходятся только Re.")

    # 2. Контроль eps_p = eps_h
    print("\n2) Контроль: частицы неотличимы от матрицы (eps_p = eps_h), обязано быть eps_h = %.4f" % EPS_H.real)
    for name, fn in MODELS.items():
        v = fn(EPS_H, EPS_H, C)
        print("   %-7s -> %8.4f%+8.4fi   ошибка %+.4f" % (name, v.real, v.imag, v.real - EPS_H.real))

    # 3. Максимумы Im eps_eff
    print("\n3) Максимум Im eps_eff (плазмонный резонанс композита)")
    peaks = {}
    for name in MODELS:
        i = max(range(len(lams)), key=lambda j: vals[name][j].imag)
        peaks[name] = (lams[i], vals[name][i].imag, vals[name][i].real)
        print("   %-7s  lambda_max = %6.1f нм   Im eps = %7.4f   Re eps = %8.4f" % (name, *peaks[name]))

    # 4. Попарная близость
    print("\n4) Попарное расхождение max|eps_A - eps_B| и max|Im eps_A - Im eps_B| на 400-750 нм")
    names = list(MODELS)
    rows = []
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            dfull = max(abs(x - y) for x, y in zip(vals[a], vals[b]))
            dim = max(abs(x.imag - y.imag) for x, y in zip(vals[a], vals[b]))
            rows.append((dfull, dim, a, b))
    for dfull, dim, a, b in sorted(rows):
        print("   %-7s vs %-7s   max|deps| = %8.4f   max|dIm| = %8.4f" % (a, b, dfull, dim))

    # 5. Значения в нескольких точках
    print("\n5) Контрольные точки")
    for lam in (450.0, 535.0, 600.0, 650.0, 700.0, 750.0):
        j = lams.index(lam)
        print("   lambda = %.0f нм, eps_Au = %8.4f%+8.4fi" % (lam, eps_p[lam].real, eps_p[lam].imag))
        for name in names:
            v = vals[name][j]
            print("      %-7s %9.4f%+9.4fi" % (name, v.real, v.imag))


if __name__ == "__main__":
    main()
