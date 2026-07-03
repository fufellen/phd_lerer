"""Валидация rcwa_2d.py (полновекторная FMM для crossed grating) против
независимых пределов - ДО какого-либо сравнения со статьями Лерера.

Тесты:
1. J1 (локальная аппроксимация A&S) против ряда Тейлора.
2. Фурье-коэффициенты диска: DC-член = доля площади.
3. Однородный предел (диаметр=0) против уже провалидированного tmm.py.
4. 1D-предел (профиль-полоска без зависимости от y) против уже
   провалидированного rcwa_1d_te.py - сильный тест всей векторной машины.
5. Сохранение энергии sum(R)+sum(T)=1 для беспотерьной решётки цилиндров.
6. Симметрия x-pol == y-pol для цилиндров на квадратной решётке при
   нормальном падении.
7. Сходимость R_00 по числу гармоник N.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from rcwa_2d import (  # noqa: E402
    bessel_j1, fourier_coeffs_cylinder, fourier_coeffs_stripe_2d,
    convolution_matrix_2d, solve_rcwa_2d,
)
from rcwa_1d_te import fourier_coeffs_strip, toeplitz_from_coeffs, solve_rcwa_te  # noqa: E402
from tmm import Layer, power_rt  # noqa: E402


def test_bessel_j1():
    """J1 по ряду Тейлора: J1(x) = sum_k (-1)^k (x/2)^{2k+1} / (k! (k+1)!)."""
    for x in (0.1, 0.5, 1.0, 2.0, 5.0, 7.9, 8.1, 12.0, 20.0):
        s, term = 0.0, x / 2.0
        k = 0
        while abs(term) > 1e-18 and k < 200:
            s += term
            k += 1
            term *= -(x / 2.0) ** 2 / (k * (k + 1))
        assert abs(bessel_j1(x) - s) < 5e-8, (x, bessel_j1(x), s)
        assert abs(bessel_j1(-x) + s) < 5e-8  # нечётность
    print("[OK] bessel_j1 совпадает с рядом Тейлора (|err|<5e-8), нечётная.")


def test_fourier_disk_dc():
    eps_incl, eps_bg = 4.0 + 0j, 1.0 + 0j
    d_over_p = 0.8
    a = fourier_coeffs_cylinder(eps_incl, eps_bg, d_over_p, N=4)
    f_area = math.pi * (d_over_p / 2) ** 2
    expected = eps_bg + (eps_incl - eps_bg) * f_area
    assert abs(a[4, 4] - expected) < 1e-12
    # симметрия диска: a[m,n]=a[n,m]=a[-m,n]
    assert abs(a[4 + 1, 4 + 2] - a[4 + 2, 4 + 1]) < 1e-15
    assert abs(a[4 + 1, 4 + 2] - a[4 - 1, 4 + 2]) < 1e-15
    print(f"[OK] Фурье диска: DC = eps_bg + delta*pi*R^2 = {expected:.4f}; симметрии выполняются.")


def test_homogeneous_limit_vs_tmm():
    """Слой 'решётки' с diameter=0 -> однородный eps_bg: R,T должны совпасть
    с TMM, вся мощность в порядке (0,0)."""
    N = 2
    lam, period, d = 600.0, 300.0, 80.0
    eps_bg = 2.1 + 0.05j
    eps_in, eps_out = 1.45 ** 2 + 0j, 1.0 + 0j

    a = fourier_coeffs_cylinder(eps_bg, eps_bg, 0.5, N)  # eps_incl=eps_bg: диск невидим
    E = convolution_matrix_2d(a, N)
    R_mn, T_mn = solve_rcwa_2d([(E, d)], eps_in, eps_out, N, period, lam, pol="y")

    import cmath
    n_bg = cmath.sqrt(eps_bg)
    if n_bg.imag < 0:
        n_bg = -n_bg
    R_tmm, T_tmm = power_rt([Layer(1.45, 0.0), Layer(n_bg, d), Layer(1.0, 0.0)], lam, 0.0, "s")

    R00, T00 = R_mn[N, N], T_mn[N, N]
    R_other = R_mn.sum() - R00
    T_other = T_mn.sum() - T00
    print(f"2D RCWA (0,0): R={R00:.6f} T={T00:.6f}   TMM: R={R_tmm:.6f} T={T_tmm:.6f}")
    assert abs(R00 - R_tmm) < 1e-6, (R00, R_tmm)
    assert abs(T00 - T_tmm) < 1e-6, (T00, T_tmm)
    assert abs(R_other) < 1e-9 and abs(T_other) < 1e-9
    print("[OK] Однородный предел: 2D RCWA == TMM, вся мощность в (0,0).")


def test_stripe_limit_vs_rcwa_1d():
    """Профиль-полоска (нет зависимости от y) при Ey-поляризации: полный
    векторный 2D-код должен воспроизвести уже провалидированный 1D TE-код."""
    N = 3
    lam, period, d = 600.0, 400.0, 150.0
    eps_incl, eps_bg = 4.0 + 0.2j, 1.0 + 0j
    eps_in, eps_out = 1.5 ** 2 + 1e-9j, 1.0 + 1e-9j
    fill = 0.4

    a2 = fourier_coeffs_stripe_2d(eps_incl, eps_bg, fill, N)
    E2 = convolution_matrix_2d(a2, N)
    R_mn, T_mn = solve_rcwa_2d([(E2, d)], eps_in, eps_out, N, period, lam, pol="y")

    c1 = fourier_coeffs_strip(eps_incl, eps_bg, fill, N)
    T1 = toeplitz_from_coeffs(c1)
    R1_harm, T1_harm = solve_rcwa_te([(T1, d)], eps_in, eps_out, N, period, lam, theta_i_deg=0.0)

    # в 2D-коде порядки (m, n): при профиле без y-зависимости всё в n=0
    R_m = R_mn[:, N]
    T_m = T_mn[:, N]
    off_axis = R_mn.sum() - R_m.sum() + T_mn.sum() - T_m.sum()
    max_dR = np.max(np.abs(R_m - R1_harm))
    max_dT = np.max(np.abs(T_m - T1_harm))
    print(f"1D-предел: max|dR_m|={max_dR:.2e}, max|dT_m|={max_dT:.2e}, "
          f"мощность вне n=0: {off_axis:.2e}")
    assert max_dR < 1e-8 and max_dT < 1e-8
    assert abs(off_axis) < 1e-10
    print("[OK] 1D-предел: 2D-код повторяет rcwa_1d_te по каждому порядку m.")


def test_energy_conservation_lossless_cylinders():
    N = 3
    period, d = 500.0, 100.0
    eps_incl, eps_bg = 4.0 + 0j, 1.0 + 0j  # без потерь
    eps_in, eps_out = 1.0 + 1e-9j, 1.45 ** 2 + 1e-9j
    a = fourier_coeffs_cylinder(eps_incl, eps_bg, 0.8, N)
    E = convolution_matrix_2d(a, N)
    max_imbalance = 0.0
    for lam in np.arange(450.0, 900.0, 37.0):
        R_mn, T_mn = solve_rcwa_2d([(E, d)], eps_in, eps_out, N, period, float(lam), pol="y")
        max_imbalance = max(max_imbalance, abs(R_mn.sum() + T_mn.sum() - 1.0))
    print(f"[energy] max |sum(R)+sum(T)-1| (беспотерьные цилиндры): {max_imbalance:.2e}")
    assert max_imbalance < 1e-6
    print("[OK] Энергия сохраняется для беспотерьной решётки цилиндров.")


def test_xy_polarization_symmetry():
    """Цилиндры на квадратной решётке при нормальном падении: x- и y-поляризация
    эквивалентны (поворот на 90 градусов - симметрия структуры)."""
    N = 3
    lam, period, d = 600.0, 500.0, 100.0
    eps_incl, eps_bg = 3.0 + 1.5j, 1.0 + 0j
    eps_in, eps_out = 1.0 + 1e-9j, 1.45 ** 2 + 1e-9j
    a = fourier_coeffs_cylinder(eps_incl, eps_bg, 0.8, N)
    E = convolution_matrix_2d(a, N)
    R_x, T_x = solve_rcwa_2d([(E, d)], eps_in, eps_out, N, period, lam, pol="x")
    R_y, T_y = solve_rcwa_2d([(E, d)], eps_in, eps_out, N, period, lam, pol="y")
    dR = abs(R_x.sum() - R_y.sum())
    dT = abs(T_x.sum() - T_y.sum())
    # порядки при повороте меняются местами (m,n)->(n,m)
    d_orders = np.max(np.abs(R_x - R_y.T)) + np.max(np.abs(T_x - T_y.T))
    print(f"[symmetry] |R_x-R_y|={dR:.2e}, |T_x-T_y|={dT:.2e}, порядки (транспонирование): {d_orders:.2e}")
    assert dR < 1e-9 and dT < 1e-9 and d_orders < 1e-9
    print("[OK] x-pol == y-pol для цилиндров (суммы и порядки с точностью до транспонирования).")


def test_convergence_with_N():
    lam, period, d = 600.0, 500.0, 100.0
    eps_incl, eps_bg = 3.0 + 1.5j, 1.0 + 0j
    eps_in, eps_out = 1.0 + 1e-9j, 1.45 ** 2 + 1e-9j
    prev = None
    for N in (2, 3, 4, 5):
        a = fourier_coeffs_cylinder(eps_incl, eps_bg, 0.8, N)
        E = convolution_matrix_2d(a, N)
        R_mn, T_mn = solve_rcwa_2d([(E, d)], eps_in, eps_out, N, period, lam, pol="y")
        R, T = R_mn.sum(), T_mn.sum()
        P = 1.0 - R - T
        print(f"[convergence] N={N} (гармоник {(2*N+1)**2}): R={R:.5f} T={T:.5f} P={P:.5f}")
        if prev is not None:
            assert abs(P - prev) < 0.05, "P должен сходиться при увеличении N"
        prev = P
    print("[OK] P сходится с ростом N (правило Лорана, умеренный контраст).")


def main():
    test_bessel_j1()
    test_fourier_disk_dc()
    test_homogeneous_limit_vs_tmm()
    test_stripe_limit_vs_rcwa_1d()
    test_energy_conservation_lossless_cylinders()
    test_xy_polarization_symmetry()
    test_convergence_with_N()
    print("\nВсе проверки RCWA (2D, crossed grating) пройдены.\n")


if __name__ == "__main__":
    main()
