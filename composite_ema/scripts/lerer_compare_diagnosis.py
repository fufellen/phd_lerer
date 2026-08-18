"""Диагностика сравнения трёх моделей EMA на данных самого Лерера (luxpop).

Повод: письмо 2026-08-18 «Сравнение 3 формул расчета эпс композита» - Бруггеман
и Максвелл-Гарнетт с поправкой MLWA дают близкие спектры потерь решётки, а
Клаузиус-Моссотти - сильно отличающийся; вопрос, что не так.

Скрипт считает eps_eff по трём моделям на двух наборах оптических постоянных:
 - luxpop (таблицы Au_.c / Ag_.c из базы Лерера, шаг 1 нм);
 - Johnson & Christy 1972 (используются в докладе LFPM-2026).

Параметры взяты из его же входного файла: n_m = 1,77, C = 0,10,
R_p = 30 (как передано в вызове), плюс контрольный вариант R_p = 15
на случай, если 30 нм имелись в виду как диаметр.
"""

from __future__ import annotations

import cmath
import math
import re
from pathlib import Path

import composite_ema as ce

LUXPOP = Path(__file__).resolve().parent.parent / "source_c" / "luxpop_lerer"
ROW = re.compile(r"Lam\[(\d+)\]\s*=\s*([\d.]+);\s*n\[\1\]\s*=\s*([\d.]+);\s*k\[\1\]\s*=\s*([\d.]+);")


def load_luxpop(filename: str) -> dict[int, complex]:
    """Прочитать таблицу n,k из C-исходника и перевести в конвенцию Лерера."""
    text = (LUXPOP / filename).read_text(encoding="latin-1")
    table: dict[int, complex] = {}
    for m in ROW.finditer(text):
        idx = int(m.group(1))
        n, k = float(m.group(3)), float(m.group(4))
        table[idx] = complex(n * n - k * k, -2.0 * n * k)
    return table


def eps_luxpop(table: dict[int, complex], lam: float) -> complex:
    """Воспроизводит выборку значения так, как это делает код Лерера."""
    m = int(lam + 1.0e-7)
    return table[m]


# --- три модели, конвенция Лерера exp(+i w t) ------------------------------


def mg(eps_m: complex, eps_p: complex, C: float) -> complex:
    K = (eps_p - eps_m) / (eps_p + 2.0 * eps_m)
    return eps_m * (1.0 + 2.0 * C * K) / (1.0 - C * K)


def bruggeman(eps_m: complex, eps_p: complex, C: float) -> complex:
    b = (3.0 * C - 1.0) * eps_p + (2.0 - 3.0 * C) * eps_m
    disc = cmath.sqrt(b * b + 8.0 * eps_p * eps_m)
    r1, r2 = (b + disc) / 4.0, (b - disc) / 4.0
    return r1 if abs(r1 - eps_m) <= abs(r2 - eps_m) else r2


def mlwa(eps_m: complex, eps_p: complex, C: float, lam: float, R: float) -> complex:
    eh, ep = eps_m.conjugate(), eps_p.conjugate()
    n_h = cmath.sqrt(eh).real
    x = 2.0 * math.pi * n_h * R / lam
    K = (ep - eh) / (ep + 2.0 * eh)
    S = C * K / (1.0 - x * x * K - 1j * (2.0 / 3.0) * x**3 * K)
    return (eh * (1.0 + 2.0 * S) / (1.0 - S)).conjugate()


def peak(fn, lo: float = 400.0, hi: float = 751.0) -> tuple[float, float]:
    best_lam, best = lo, 0.0
    lam = lo
    while lam <= hi:
        v = abs(fn(lam).imag)
        if v > best:
            best, best_lam = v, lam
        lam += 1.0
    return best_lam, best


def frohlich(eps_fn, eps_m_re: float, lo: float = 400.0, hi: float = 751.0):
    prev_lam, prev = None, None
    lam = lo
    while lam <= hi:
        val = eps_fn(lam).real + 2.0 * eps_m_re
        if prev is not None and prev * val < 0:
            return prev_lam + (lam - prev_lam) * abs(prev) / (abs(prev) + abs(val))
        prev, prev_lam = val, lam
        lam += 1.0
    return None


def main() -> None:
    n_m = 1.77
    eps_m = complex(n_m**2, 0.0)
    C = 0.10

    au_lux = load_luxpop("Au_.c")
    ag_lux = load_luxpop("Ag_.c")
    print("Таблицы luxpop прочитаны: Au %d узлов, Ag %d узлов"
          % (len(au_lux), len(ag_lux)))

    sources = {
        "luxpop (Лерер)": {
            "Au": lambda l: eps_luxpop(au_lux, l),
            "Ag": lambda l: eps_luxpop(ag_lux, l),
        },
        "Johnson-Christy": {
            "Au": ce.eps_au,
            "Ag": ce.eps_ag,
        },
    }

    for src_name, mats in sources.items():
        print()
        print("=" * 78)
        print("Источник оптических постоянных: %s" % src_name)
        print("=" * 78)
        for mat, eps_fn in mats.items():
            lam_f = frohlich(eps_fn, eps_m.real)
            print("\n  %s, условие Фрёлиха Re eps_p + 2 eps_m = 0: %s"
                  % (mat, ("%.0f нм" % lam_f) if lam_f else "нет в 400-751 нм"))
            rows = [
                ("Клаузиус-Моссотти", lambda l, f=eps_fn: mg(eps_m, f(l), C)),
                ("Бруггеман", lambda l, f=eps_fn: bruggeman(eps_m, f(l), C)),
                ("MLWA, R_p=30 нм", lambda l, f=eps_fn: mlwa(eps_m, f(l), C, l, 30.0)),
                ("MLWA, R_p=15 нм", lambda l, f=eps_fn: mlwa(eps_m, f(l), C, l, 15.0)),
                ("MLWA, R_p=10 нм", lambda l, f=eps_fn: mlwa(eps_m, f(l), C, l, 10.0)),
            ]
            print("    %-20s %10s %12s %10s" % ("модель", "пик, нм", "|Im eps|", "x в пике"))
            for name, fn in rows:
                lam_p, val = peak(fn)
                x = 2.0 * math.pi * n_m * (30.0 if "30" in name else 15.0 if "15" in name else 10.0) / lam_p
                x_str = "%.2f" % x if "MLWA" in name else "-"
                print("    %-20s %10.0f %12.3f %10s" % (name, lam_p, val, x_str))

    # Насколько сами eps металла различаются между базами
    print()
    print("=" * 78)
    print("Расхождение оптических постоянных: luxpop против Johnson-Christy")
    print("=" * 78)
    print("  %6s | %-26s | %-26s" % ("нм", "Au: luxpop / J&C", "Ag: luxpop / J&C"))
    for lam in (400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0, 750.0):
        a1, a2 = eps_luxpop(au_lux, lam), ce.eps_au(lam)
        g1, g2 = eps_luxpop(ag_lux, lam), ce.eps_ag(lam)
        print("  %6.0f | %+7.2f%+7.2fj / %+6.2f%+6.2fj | %+7.2f%+7.2fj / %+6.2f%+6.2fj"
              % (lam, a1.real, a1.imag, a2.real, a2.imag,
                 g1.real, g1.imag, g2.real, g2.imag))

    # Проверка интерполяции в коде Лерера
    print()
    print("=" * 78)
    print("Проверка узла интерполяции в Au_()/Ag_()")
    print("=" * 78)
    print("  Код: n_ = n[m+1]*(Lam[m]-m_) + n[m]*(m_+1-Lam[m]), где Lam[m] == m.")
    print("  Значит Lam[m]-m_ == 0 и m_+1-Lam[m] == 1, то есть n_ == n[m] всегда.")
    for lam in (600.0, 600.5, 600.9):
        v = eps_luxpop(au_lux, lam)
        print("    Au, lam=%.1f -> eps = %+.4f%+.4fj (значение узла %d)"
              % (lam, v.real, v.imag, int(lam + 1e-7)))
    print("  Линейная интерполяция фактически не выполняется: берётся нижний узел.")
    print("  При шаге таблицы 1 нм эффект мал, но заявленной интерполяции нет.")


if __name__ == "__main__":
    main()
