"""Проверка формул файла ФДП_COMP против рабочего кода репозитория.

Каждая формула реализована ЗАНОВО, буквально в том виде, в каком она записана
в ФДП_COMP (конвенция Лерера exp(+i w t), eps = eps' - i eps''), и сверяется
с composite_ema.py / compare_effective_models.py.

Цель - не «ещё один расчёт», а гарантия, что текст ФДП не расходится с кодом,
который Лерер вставляет в свой проект.
"""

from __future__ import annotations

import cmath
import math

import composite_ema as ce
import compare_effective_models as cm


# --- Реализация строго "как написано в ФДП" -------------------------------


def fdp_mg(eps_h: complex, eps_p: complex, f: float) -> complex:
    """ФДП п.2: K = (eps_p - eps_h)/(eps_p + 2 eps_h); eps = eps_h (1+2fK)/(1-fK)."""
    K = (eps_p - eps_h) / (eps_p + 2.0 * eps_h)
    return eps_h * (1.0 + 2.0 * f * K) / (1.0 - f * K)


def fdp_mg_multi(eps_h: complex, populations: list[tuple[complex, float]]) -> complex:
    """ФДП п.3: eta = sum f_j K_j; eps = eps_h (1+2 eta)/(1- eta)."""
    eta = 0j
    for eps_p, f_j in populations:
        eta += f_j * (eps_p - eps_h) / (eps_p + 2.0 * eps_h)
    return eps_h * (1.0 + 2.0 * eta) / (1.0 - eta)


def fdp_bruggeman(eps_h: complex, eps_p: complex, f: float) -> complex:
    """ФДП п.4: 2 eps^2 - b eps - eps_p eps_h = 0,
    b = (3f-1) eps_p + (2-3f) eps_h, корень ближе к eps_h."""
    b = (3.0 * f - 1.0) * eps_p + (2.0 - 3.0 * f) * eps_h
    disc = cmath.sqrt(b * b + 8.0 * eps_p * eps_h)
    r1 = (b + disc) / 4.0
    r2 = (b - disc) / 4.0
    return r1 if abs(r1 - eps_h) <= abs(r2 - eps_h) else r2


def fdp_mlwa(eps_h: complex, eps_p: complex, f: float, lam_nm: float, R_nm: float) -> complex:
    """ФДП п.5, безразмерная форма:
    x = k_h R = 2 pi n_h R / lam,
    eta = f K / (1 - x^2 K - i (2/3) x^3 K),
    eps = eps_h (1 + 2 eta)/(1 - eta).
    Формула записана в конвенции exp(-i w t): сопрягаем на входе и на выходе.
    """
    eh = eps_h.conjugate()
    ep = eps_p.conjugate()

    n_h = cmath.sqrt(eh).real
    x = 2.0 * math.pi * n_h * R_nm / lam_nm

    K = (ep - eh) / (ep + 2.0 * eh)
    eta = f * K / (1.0 - x * x * K - 1j * (2.0 / 3.0) * x**3 * K)
    eps_eff = eh * (1.0 + 2.0 * eta) / (1.0 - eta)
    return eps_eff.conjugate()


# --- Сверка ----------------------------------------------------------------


def main() -> None:
    ce.self_test()

    n_h = 1.77
    eps_h = complex(n_h**2, 0.0)
    f = 0.10
    tol = 1e-10
    worst = {"mg": 0.0, "brug": 0.0, "mlwa": 0.0}

    lam = 300.0
    while lam <= 900.0 + 1e-9:
        for mat_fn in (ce.eps_au, ce.eps_cu):
            eps_p = mat_fn(lam)

            d = abs(fdp_mg(eps_h, eps_p, f) - ce.composite_from_eps(eps_h, eps_p, f))
            worst["mg"] = max(worst["mg"], d)
            assert d < tol, ("MG", lam, d)

            d = abs(fdp_bruggeman(eps_h, eps_p, f) - cm.bruggeman_from_eps(eps_h, eps_p, f))
            worst["brug"] = max(worst["brug"], d)
            assert d < 1e-9, ("Bruggeman", lam, d)

            for R in (10.0, 30.0, 50.0):
                d = abs(
                    fdp_mlwa(eps_h, eps_p, f, lam, R)
                    - cm.mlwa_mg_from_eps(eps_h, eps_p, f, lam, R)
                )
                worst["mlwa"] = max(worst["mlwa"], d)
                assert d < 1e-9, ("MLWA", lam, R, d)
        lam += 2.0

    print("[ФДП] MG      max |diff| = %.3e" % worst["mg"])
    print("[ФДП] Бруггеман max |diff| = %.3e" % worst["brug"])
    print("[ФДП] MLWA    max |diff| = %.3e" % worst["mlwa"])

    # Предельные случаи, заявленные в ФДП п.8
    eps_p = ce.eps_au(535.0)
    assert fdp_mg(eps_h, eps_p, 0.0) == eps_h, "f=0 должно вернуть eps_h точно"
    assert abs(fdp_bruggeman(eps_h, eps_p, 0.0) - eps_h) < 1e-12, "Бруггеман при f=0"
    # один сорт частиц: многосортовая формула сводится к однокомпонентной
    assert abs(fdp_mg_multi(eps_h, [(eps_p, f)]) - fdp_mg(eps_h, eps_p, f)) < 1e-12
    # согласование с портом COMPOSITE_3
    two = [(ce.eps_au(535.0), 0.06), (ce.eps_cu(535.0), 0.04)]
    assert abs(fdp_mg_multi(eps_h, two) - ce.composite_3_from_eps(eps_h, two)) < 1e-12
    # MLWA при x -> 0 переходит в квазистатический Максвелл-Гарнетт.
    # Поправка ведущего порядка ~ x^2, поэтому проверяем именно порядок
    # сходимости: при уменьшении R вдвое расхождение падает примерно вчетверо.
    d_1 = abs(fdp_mlwa(eps_h, eps_p, f, 535.0, 1.0) - fdp_mg(eps_h, eps_p, f))
    d_2 = abs(fdp_mlwa(eps_h, eps_p, f, 535.0, 0.5) - fdp_mg(eps_h, eps_p, f))
    ratio = d_1 / d_2
    assert 3.7 < ratio < 4.3, ("порядок сходимости MLWA -> MG", ratio, d_1, d_2)
    print(
        "[ФДП] предельные случаи (f=0, один сорт, COMPOSITE_3) - OK; "
        "MLWA -> MG при x->0 со скоростью x^2 (отношение %.2f при R 1->0,5 нм)" % ratio
    )

    # Пассивность в конвенции Лерера
    lam = 300.0
    while lam <= 900.0 + 1e-9:
        for mat_fn in (ce.eps_au, ce.eps_cu):
            eps_p = mat_fn(lam)
            for val in (
                fdp_mg(eps_h, eps_p, f),
                fdp_bruggeman(eps_h, eps_p, f),
                fdp_mlwa(eps_h, eps_p, f, lam, 10.0),
                fdp_mlwa(eps_h, eps_p, f, lam, 30.0),
                fdp_mlwa(eps_h, eps_p, f, lam, 50.0),
            ):
                assert val.imag <= 1e-9, ("пассивность нарушена", lam, val)
        lam += 2.0
    print("[ФДП] пассивность Im eps_eff <= 0 на всей сетке 300-900 нм - OK")

    # Контрольная точка ФДП п.10
    eps_p = ce.eps_au(535.0)
    print("\nКонтрольная точка 535 нм, Au, n_h=1.77, f=0.10:")
    print("  Максвелл-Гарнетт  %s" % _fmt(fdp_mg(eps_h, eps_p, f)))
    print("  Бруггеман         %s" % _fmt(fdp_bruggeman(eps_h, eps_p, f)))
    print("  MLWA R=30 нм      %s" % _fmt(fdp_mlwa(eps_h, eps_p, f, 535.0, 30.0)))

    # Потери (ФДП п.6) в этой же точке
    for name, val in (
        ("Максвелл-Гарнетт", fdp_mg(eps_h, eps_p, f)),
        ("Бруггеман", fdp_bruggeman(eps_h, eps_p, f)),
        ("MLWA R=30 нм", fdp_mlwa(eps_h, eps_p, f, 535.0, 30.0)),
    ):
        N = cmath.sqrt(val.conjugate())      # в конвенцию exp(-i w t)
        n_e, k_e = N.real, abs(N.imag)
        alpha_um = 4.0 * math.pi * k_e / (535.0e-3)   # 1/мкм
        print(
            "  %-17s n=%.4f  kappa=%.4f  alpha=%.4f 1/мкм  = %.2f дБ/мкм"
            % (name, n_e, k_e, alpha_um, 10.0 * math.log10(math.e) * alpha_um)
        )

    # Условие Фрёлиха
    cross = ce.find_frohlich_crossing(ce.eps_au, eps_h.real, 400.0, 700.0)
    print("\nУсловие Фрёлиха Re eps_p + 2 eps_h = 0 (Au): %.0f нм" % cross.lam_nm)
    cross_cu = ce.find_frohlich_crossing(ce.eps_cu, eps_h.real, 400.0, 700.0)
    print("Условие Фрёлиха Re eps_p + 2 eps_h = 0 (Cu): %.0f нм" % cross_cu.lam_nm)

    print("\n[ФДП] ВСЕ ФОРМУЛЫ ФАЙЛА ФДП_COMP СОВПАДАЮТ С КОДОМ РЕПОЗИТОРИЯ.")


def _fmt(z: complex) -> str:
    return "%+.4f%+.4fj" % (z.real, z.imag)


if __name__ == "__main__":
    main()
