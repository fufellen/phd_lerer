"""
Флагманская фигура нелинейной статьи: спектр поглощения P(lambda) решётки с
золотыми ПНЧ при нескольких интенсивностях I. Связывает нелинейность золота
(ветка A, феноменологическая chi^(3)) с реальной метрикой метаповерхности -
поглощением P = 1 - R - T.

Схема расчёта на каждой (lambda, I):
1. Ветка A: нелинейная проницаемость частицы eps_p(I) самосогласованно через
   локальное поле (модуль nonlinear_au_sweep.py), затем eps_eff(I) композита
   по Максвеллу-Гарнетту (composite_ema.composite_from_eps).
2. RCWA (1D-приближение решётки, rcwa_1d_te): eps_eff(I) подставляется и в
   сплошной слой с ПНЧ, и в материал цилиндров-полосок; решается дифракция,
   считается P(lambda, I).

Геометрия - как в линейном воспроизведении (reproduce_fig2_au_vs_cu.py):
подложка n=1.45, слой с ПНЧ h=100 нм и решётка полосок h=100 нм из композита
n_host=1.77, C=10%, период 500 нм, полоски fill=0.8 (диаметр 400 нм), TE,
нормальное падение.

ВАЖНЫЕ ОГОВОРКИ (статус: diagnostic / screening).
- Расцепленное (decoupled) приближение: нелинейность управляется ПАДАЮЩЕЙ
  интенсивностью I, а не самосогласованным локальным полем внутри решётки.
  То есть eps_eff(I) считается так, будто композит видит однородную I, после
  чего RCWA линейна с этой eps_eff. Полная нелинейная RCWA (поле в решётке
  меняет eps, eps меняет поле) здесь НЕ решается - это скрининг, показывающий
  порядок эффекта и направление сдвига/выцветания резонанса.
- chi^(3) взят как один benchmark (Hache 1988, ps/SPR, |chi_m|~5e-8 esu,
  фаза 80 град) - представитель пс/нс режима у резонанса; реальный разброс
  литературных chi^(3) на порядки (см. таблицу benchmarks рабочей заметки).
- 1D-приближение решётки (полоски вместо 2D-цилиндров), как в линейной работе.

Запуск (из корня репозитория):
    python grating_2d_rcwa/scripts/rcwa_nonlinear_au_absorption.py
Выход: grating_2d_rcwa/results/rcwa_nonlinear_au_absorption.{csv,png}
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "composite_ema" / "scripts"))

from rcwa_1d_te import fourier_coeffs_strip, toeplitz_from_coeffs, solve_rcwa_te  # noqa: E402
import composite_ema as ema  # noqa: E402
import nonlinear_au_sweep as nl  # noqa: E402

# --- геометрия (как в линейном воспроизведении рис. 2) ---
N_HOST = 1.77
C_FILL = 0.10
N_SUB = 1.45
PERIOD = 500.0
DIAMETER = 400.0
FILL = DIAMETER / PERIOD
H_CYL = 100.0
H_LAYER = 100.0
N_HARM = 6
EPS_IN = 1.0 + 1e-9j              # воздух сверху
EPS_OUT = N_SUB ** 2 + 1e-9j      # подложка

# benchmark chi^(3) для ветки A (Hache 1988, ps/SPR) в конвенции Лерера
BENCH = nl.BENCHMARKS[0]
EPS_H = complex(N_HOST ** 2, 0.0)

# интенсивности (Вт/см^2): 0 = линейный; остальные в пс/нс возмущательной ветке
INTENSITIES = [0.0, 1.0e6, 5.0e6, 2.0e7]


def eps_eff_nonlinear(lam_nm: float, i_wcm2: float) -> tuple[complex, str]:
    """eps_eff(lambda, I) композита с Au ПНЧ (конвенция composite_ema eps'-i eps'')."""
    eps_p, conv = nl.eps_p_nonlinear(lam_nm, i_wcm2, BENCH.chi_lerer_si, EPS_H)
    status = nl.perturbation_status(eps_p, lam_nm) if conv else "no_convergence"
    return ema.composite_from_eps(EPS_H, eps_p, C_FILL), status


def rcwa_P(lam_nm: float, eps_eff_leconv: complex) -> tuple[float, float, float]:
    """P = 1 - R - T решётки при данной eps_eff (в конвенции Лерера)."""
    eps_std = eps_eff_leconv.conjugate()  # -> стандартная фотонная конвенция
    coeffs = fourier_coeffs_strip(eps_std, 1.0 + 0j, FILL, N_HARM)
    T_cyl = toeplitz_from_coeffs(coeffs)
    coeffs_flat = fourier_coeffs_strip(eps_std, eps_std, 0.0, N_HARM)
    T_flat = toeplitz_from_coeffs(coeffs_flat)
    layers = [(T_cyl, H_CYL), (T_flat, H_LAYER)]
    R_harm, T_harm = solve_rcwa_te(layers, EPS_IN, EPS_OUT, N_HARM, PERIOD, lam_nm, theta_i_deg=0.0)
    R = float(np.sum(R_harm).real)
    T = float(np.sum(T_harm).real)
    return R, T, 1.0 - R - T


def self_test() -> None:
    # I=0 нелинейной ветки должно совпасть с линейной eps_eff (и, значит, P).
    lam = 535.0
    eps_lin = ema.composite_from_eps(EPS_H, ema.eps_au(lam), C_FILL)
    eps_nl0, _ = eps_eff_nonlinear(lam, 0.0)
    assert abs(eps_nl0 - eps_lin) < 1e-9, (eps_nl0, eps_lin)
    P_lin = rcwa_P(lam, eps_lin)[2]
    P_nl0 = rcwa_P(lam, eps_nl0)[2]
    assert abs(P_lin - P_nl0) < 1e-9
    print(f"[self_test] OK: I=0 воспроизводит линейный расчёт (P={P_lin:.4f} на {lam:.0f} нм).")


def main() -> None:
    self_test()
    lo, hi = ema.MATERIAL_RANGE_NM[4]  # Au
    lam_grid = np.arange(max(400.0, lo), min(900.0, hi), 5.0)

    curves = {}   # I -> (P array)
    summary = []
    rows = []
    for i_wcm2 in INTENSITIES:
        P_arr, worst = [], "ok"
        for lam in lam_grid:
            eps_eff, status = eps_eff_nonlinear(float(lam), i_wcm2)
            R, T, P = rcwa_P(float(lam), eps_eff)
            P_arr.append(P)
            if status != "ok":
                worst = status
            rows.append((i_wcm2, float(lam), R, T, P, status))
        P_arr = np.array(P_arr)
        assert np.all(P_arr > -1e-6), f"P<0 при I={i_wcm2}"
        pk = int(np.argmax(P_arr))
        curves[i_wcm2] = P_arr
        summary.append((i_wcm2, lam_grid[pk], P_arr[pk], worst))
        print(f"[I={i_wcm2:.0e}] пик P={P_arr[pk]:.3f} на {lam_grid[pk]:.0f} нм [ветка A: {worst}]")

    # сдвиг и изменение пика относительно линейного
    P0 = curves[0.0]
    pk0 = int(np.argmax(P0))
    print(f"[итог] линейный пик: P={P0[pk0]:.3f} на {lam_grid[pk0]:.0f} нм; "
          f"при I={INTENSITIES[-1]:.0e}: dP_пик={curves[INTENSITIES[-1]][pk0]-P0[pk0]:+.3f}")

    out_dir = HERE.parent / "results"
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / "rcwa_nonlinear_au_absorption.csv"
    with csv_path.open("w", encoding="utf-8") as fh:
        fh.write("I_wcm2,lambda_nm,R,T,P,branchA_status\n")
        for i_wcm2, lam, R, T, P, st in rows:
            fh.write(f"{i_wcm2},{lam},{R},{T},{P},{st}\n")
    print(f"CSV записан: {csv_path}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax, axd) = plt.subplots(2, 1, figsize=(8.0, 7.6), sharex=True)
    colors = plt.cm.viridis(np.linspace(0.0, 0.85, len(INTENSITIES)))
    for color, i_wcm2 in zip(colors, INTENSITIES):
        label = "линейный (I=0)" if i_wcm2 == 0 else f"I = {i_wcm2:.0e} Вт/см²"
        ax.plot(lam_grid, curves[i_wcm2], color=color, label=label)
        if i_wcm2 != 0.0:
            axd.plot(lam_grid, curves[i_wcm2] - P0, color=color, label=label)
    ax.set_ylabel("P = 1 - R - T (поглощение)")
    ax.set_title("RCWA: поглощение решётки с Au ПНЧ при росте интенсивности\n"
                 "(ветка A, Hache 1988; 1D-решётка; decoupled screening)", fontsize=10)
    ax.legend(); ax.grid(True, alpha=0.25)
    axd.axhline(0, color="gray", linewidth=0.6)
    axd.set_xlabel("длина волны, нм")
    axd.set_ylabel("ΔP = P(I) − P(0)")
    axd.set_title("Изменение поглощения относительно линейного (зум эффекта)", fontsize=9)
    axd.legend(fontsize=8); axd.grid(True, alpha=0.25)
    fig.tight_layout()
    png_path = out_dir / "rcwa_nonlinear_au_absorption.png"
    fig.savefig(png_path, dpi=150)
    print(f"PNG записан: {png_path}")


if __name__ == "__main__":
    main()
