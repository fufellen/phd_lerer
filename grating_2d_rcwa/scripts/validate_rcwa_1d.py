"""Валидация rcwa_1d_te.py: пределы и энергетический баланс, независимо
от статей Лерера (см. rcwa_1d_te.py про метод и его происхождение)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "composite_ema" / "scripts"))

from rcwa_1d_te import (  # noqa: E402
    fourier_coeffs_strip, toeplitz_from_coeffs, layer_eigenmodes_te,
    homogeneous_kz, solve_rcwa_te,
)
from tmm import Layer, power_rt  # noqa: E402
import composite_ema as ema  # noqa: E402


def test_fourier_dc_term():
    eps_incl, eps_bg, fill = 4.0 + 0j, 2.0 + 0j, 0.3
    c = fourier_coeffs_strip(eps_incl, eps_bg, fill, N=5)
    dc = c[5]  # N=5 -> индекс 5 = гармоника 0
    expected = fill * eps_incl + (1 - fill) * eps_bg
    assert abs(dc - expected) < 1e-12, (dc, expected)
    print(f"[OK] Fourier DC-член = {dc:.4f} = f*eps_incl+(1-f)*eps_bg = {expected:.4f}")


def test_homogeneous_eigenmodes_match_direct():
    N = 4
    eps = 2.25 + 0.1j
    coeffs = np.zeros(2 * N + 1, dtype=complex)
    coeffs[N] = eps  # только DC член - однородная среда, без решётки
    T = toeplitz_from_coeffs(coeffs)
    kx0 = 0.3
    Gx_over_k0 = 0.8  # безразмерно: (2*pi/period)/k0, произвольное для теста значение
    harmonics = np.arange(-N, N + 1)
    Kx = np.diag(kx0 + harmonics * Gx_over_k0)
    modes = layer_eigenmodes_te(T, Kx)
    kz_direct = homogeneous_kz(eps, Kx)
    kz_from_modes = sorted(np.diag(modes.Q), key=lambda v: (v.real, v.imag))
    kz_direct_sorted = sorted(kz_direct, key=lambda v: (v.real, v.imag))
    for a, b in zip(kz_from_modes, kz_direct_sorted):
        assert abs(a - b) < 1e-9, (a, b, kz_from_modes, kz_direct_sorted)
    print("[OK] Собственные моды однородного слоя совпадают с прямым расчётом kz по гармоникам.")


def test_rcwa_matches_tmm_for_uniform_layer():
    """Решётка с fill=0 (нет неоднородности) должна давать те же R,T, что и TMM
    для соответствующего однородного слоя, при этом ВСЯ мощность должна
    оставаться в (0,0)-гармонике (никакой дифракции без реальной решётки)."""
    N = 3
    lam = 600.0
    period = 300.0
    d = 80.0
    eps_bg = 2.1 + 0.05j
    eps_in, eps_out = 1.45 ** 2 + 0j, 1.0 + 0j

    coeffs = fourier_coeffs_strip(eps_bg, eps_bg, 0.0, N)  # fill=0 -> eps(x)=eps_bg everywhere
    T = toeplitz_from_coeffs(coeffs)
    R_harm, T_harm = solve_rcwa_te([(T, d)], eps_in, eps_out, N, period, lam, theta_i_deg=0.0)

    layers_tmm = [Layer(1.45, 0.0), Layer(complex(eps_bg) ** 0.5, d), Layer(1.0, 0.0)]
    # branch fix: use consistent sqrt with Im>=0
    import cmath
    n_bg = cmath.sqrt(eps_bg)
    if n_bg.imag < 0:
        n_bg = -n_bg
    layers_tmm = [Layer(1.45, 0.0), Layer(n_bg, d), Layer(1.0, 0.0)]
    R_tmm, T_tmm = power_rt(layers_tmm, lam, 0.0, "s")

    R0 = R_harm[N]
    T0 = T_harm[N]
    R_other = np.sum(R_harm) - R0
    T_other = np.sum(T_harm) - T0
    print(f"RCWA (0,0): R={R0:.6f} T={T0:.6f}   TMM: R={R_tmm:.6f} T={T_tmm:.6f}")
    print(f"RCWA прочие гармоники: R_other={R_other:.2e} T_other={T_other:.2e} (должны быть ~0)")
    assert abs(R0 - R_tmm) < 1e-6, (R0, R_tmm)
    assert abs(T0 - T_tmm) < 1e-6, (T0, T_tmm)
    assert abs(R_other) < 1e-9 and abs(T_other) < 1e-9
    print("[OK] RCWA с fill=0 (нет решётки) точно совпадает с TMM для однородного слоя; вся мощность в (0,0).")


def test_energy_conservation_lossless_grating():
    N = 5
    lam = 600.0
    period = 400.0  # достаточно большой период, чтобы были пропагирующие высшие порядки
    d = 150.0
    eps_incl, eps_bg = 4.0 + 0j, 1.0 + 0j  # без потерь
    eps_in, eps_out = 1.5 ** 2 + 0j, 1.0 + 0j
    fill = 0.4

    coeffs = fourier_coeffs_strip(eps_incl, eps_bg, fill, N)
    T = toeplitz_from_coeffs(coeffs)
    max_imbalance = 0.0
    for lam_scan in np.arange(450.0, 900.0, 23.0):
        R_harm, T_harm = solve_rcwa_te([(T, d)], eps_in, eps_out, N, period, float(lam_scan), theta_i_deg=0.0)
        total = R_harm.sum() + T_harm.sum()
        max_imbalance = max(max_imbalance, abs(total.real - 1.0))
    print(f"[energy] максимальный |sum(R)+sum(T)-1| по сетке длин волн (беспотерьная решётка): {max_imbalance:.2e}")
    assert max_imbalance < 1e-6, max_imbalance
    print("[OK] Энергия сохраняется (R+T=1 суммарно по всем гармоникам) для беспотерьной решётки.")


def test_convergence_with_N():
    lam = 600.0
    period = 400.0
    d = 150.0
    eps_incl, eps_bg = 4.0 + 0j, 1.0 + 0j
    eps_in, eps_out = 1.5 ** 2 + 0j, 1.0 + 0j
    fill = 0.4
    prev_R0 = None
    for N in (2, 4, 6, 9):
        coeffs = fourier_coeffs_strip(eps_incl, eps_bg, fill, N)
        T = toeplitz_from_coeffs(coeffs)
        R_harm, T_harm = solve_rcwa_te([(T, d)], eps_in, eps_out, N, period, lam, theta_i_deg=0.0)
        R0 = R_harm[N]
        print(f"[convergence] N={N}: R_0 = {R0:.6f}  (sum R+T = {R_harm.sum()+T_harm.sum():.6f})")
        if prev_R0 is not None:
            assert abs(R0 - prev_R0) < 0.05, "R_0 должен сходиться при увеличении N"
        prev_R0 = R0
    print("[OK] R_0 сходится с ростом числа гармоник N.")


def main():
    test_fourier_dc_term()
    test_homogeneous_eigenmodes_match_direct()
    test_rcwa_matches_tmm_for_uniform_layer()
    test_energy_conservation_lossless_grating()
    test_convergence_with_N()
    print("\nВсе проверки RCWA (1D, TE) пройдены.\n")


if __name__ == "__main__":
    main()
