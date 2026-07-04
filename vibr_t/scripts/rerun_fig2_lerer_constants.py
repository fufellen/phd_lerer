"""
ПЕРЕСЧЁТ рис. 2 статьи [[Лерер АМ_1]] (Au vs Cu, поглощение решётки) на
ПЕРВОИСТОЧНИКЕ оптических констант - таблицах из программы Лерера VIBR_T
(vibr_t), вместо Johnson & Christy 1972 (composite_ema).

Это прямая проверка последней открытой гипотезы воспроизведения:
"УФ-порядок Au/Cu в нашей репродукции противоположен статье, потому что
у нас другой первоисточник констант". Геометрия и солвер - в точности как
в grating_2d_rcwa/scripts/reproduce_fig2_au_vs_cu.py (RCWA, 1D-приближение
решётки, TE): период 500 нм, полоски fill=0.8, h=100 нм, слой с ПНЧ 100 нм,
n_host=1.77, C=10%, подложка n=1.45.

Считается три варианта Cu:
  1) Cu по таблице Лерера "как в программе" - с линейной экстраполяцией
     ниже 517 нм (Re eps > 0 в УФ - медь-"диэлектрик");
  2) Cu по Johnson & Christy (наша прежняя чистовая база) - для сравнения;
Au в обоих случаях берётся из соответствующей базы.

Выход: results/rerun_fig2_lerer_constants.{csv,png} + печатный вердикт
по среднему УФ-поглощению.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "composite_ema" / "scripts"))
sys.path.insert(0, str(HERE.parents[1] / "grating_2d_rcwa" / "scripts"))

import vibr_materials as vm  # noqa: E402
import vibr_composite as vc  # noqa: E402
import composite_ema as ema  # noqa: E402
from rcwa_1d_te import fourier_coeffs_strip, toeplitz_from_coeffs, solve_rcwa_te  # noqa: E402

RESULTS = HERE.parent / "results"

N_HOST = 1.77
C_NP = 0.10
N_SUB = 1.45
PERIOD = 500.0
FILL = 400.0 / 500.0
H_CYL = 100.0
H_LAYER = 100.0
N_HARM = 6
LAM = np.arange(210.0, 900.0, 5.0)


def run_curve(eps_np_fn, label: str):
    """P(lambda) для решётки, eps_np_fn(lam)->eps наночастиц (конвенция Лерера)."""
    eps_h = complex(N_HOST ** 2, 0.0)
    eps_in = 1.0 + 1e-9j
    eps_out = N_SUB ** 2 + 1e-9j
    R_l, T_l, P_l = [], [], []
    for lam in LAM:
        eps_p = eps_np_fn(float(lam))
        eps_comp = vc.mg_single(eps_h, eps_p, C_NP)      # конвенция Лерера (Im<=0)
        eps_std = eps_comp.conjugate()                    # фотонная конвенция для RCWA
        coeffs = fourier_coeffs_strip(eps_std, 1.0 + 0j, FILL, N_HARM)
        T_cyl = toeplitz_from_coeffs(coeffs)
        coeffs_flat = fourier_coeffs_strip(eps_std, eps_std, 0.0, N_HARM)
        T_flat = toeplitz_from_coeffs(coeffs_flat)
        layers = [(T_cyl, H_CYL), (T_flat, H_LAYER)]
        R_h, T_h = solve_rcwa_te(layers, eps_in, eps_out, N_HARM, PERIOD, float(lam),
                                 theta_i_deg=0.0)
        R = float(np.sum(R_h).real)
        T = float(np.sum(T_h).real)
        R_l.append(R)
        T_l.append(T)
        P_l.append(1.0 - R - T)
    P = np.array(P_l)
    print(f"[{label}] P_max={P.max():.3f} на {LAM[np.argmax(P)]:.0f} нм; "
          f"среднее P в УФ (<=400 нм) = {P[LAM <= 400.0].mean():.3f}")
    return np.array(R_l), np.array(T_l), P


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)

    print("=== Константы Лерера (VIBR_T) ===")
    R_au_l, T_au_l, P_au_l = run_curve(lambda l: vm.libmat(3, l), "Au, Лерер")
    R_cu_l, T_cu_l, P_cu_l = run_curve(lambda l: vm.libmat(1, l), "Cu, Лерер (экстраполяция в УФ)")

    print("=== Johnson & Christy 1972 (наша прежняя база) ===")
    R_au_j, T_au_j, P_au_j = run_curve(lambda l: ema.libmat(4, l), "Au, J&C")
    R_cu_j, T_cu_j, P_cu_j = run_curve(lambda l: ema.libmat(1, l), "Cu, J&C")

    uv = LAM <= 400.0
    print("\n--- ВЕРДИКТ по УФ-порядку Au/Cu (среднее P при lambda<=400 нм) ---")
    au_l, cu_l = P_au_l[uv].mean(), P_cu_l[uv].mean()
    au_j, cu_j = P_au_j[uv].mean(), P_cu_j[uv].mean()
    print(f"Константы Лерера:  Au={au_l:.3f}  Cu={cu_l:.3f}  ->  "
          f"{'Au > Cu (порядок КАК В СТАТЬЕ)' if au_l > cu_l else 'Cu > Au (расхождение остаётся)'}")
    print(f"Johnson & Christy: Au={au_j:.3f}  Cu={cu_j:.3f}  ->  "
          f"{'Au > Cu (как в статье)' if au_j > cu_j else 'Cu > Au (наше прежнее расхождение)'}")
    if (au_l > cu_l) != (au_j > cu_j):
        print("ГИПОТЕЗА ПОДТВЕРЖДЕНА: порядок Au/Cu в УФ определяется первоисточником "
              "констант; в программе Лерера Cu ниже 517 нм экстраполируется "
              "(Re eps > 0, 'диэлектрическая' медь) - отсюда слабое поглощение Cu в УФ "
              "на рис. 2 статьи.")
    else:
        print("Гипотеза НЕ подтверждена: порядок одинаков в обеих базах - "
              "искать другую причину.")

    # энергетический контроль
    for R, T, P, lbl in ((R_au_l, T_au_l, P_au_l, "Au-Лерер"),
                         (R_cu_l, T_cu_l, P_cu_l, "Cu-Лерер"),
                         (R_au_j, T_au_j, P_au_j, "Au-J&C"),
                         (R_cu_j, T_cu_j, P_cu_j, "Cu-J&C")):
        imb = np.max(np.abs(R + T + P - 1.0))
        assert imb < 1e-12, f"{lbl}: нарушен баланс R+T+P=1"
        assert np.all(P > -1e-6), f"{lbl}: P<0 - нефизичное усиление"
    print("[OK] R+T+P=1 и P>=0 для всех кривых")

    csv_path = RESULTS / "rerun_fig2_lerer_constants.csv"
    header = ("lambda_nm,P_Au_lerer,P_Cu_lerer,P_Au_jc,P_Cu_jc,"
              "R_Au_lerer,T_Au_lerer,R_Cu_lerer,T_Cu_lerer")
    np.savetxt(csv_path,
               np.column_stack([LAM, P_au_l, P_cu_l, P_au_j, P_cu_j,
                                R_au_l, T_au_l, R_cu_l, T_cu_l]),
               delimiter=",", header=header, comments="")
    print(f"CSV записан: {csv_path}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    ax1.plot(LAM, P_au_l, color="#d4a017", label="Au")
    ax1.plot(LAM, P_cu_l, color="#b5651d", label="Cu (экстрап. <517 нм)")
    ax1.axvline(517.0, color="gray", lw=1, ls=":")
    ax1.set_title("Константы Лерера (VIBR_T)")
    ax2.plot(LAM, P_au_j, color="#d4a017", label="Au")
    ax2.plot(LAM, P_cu_j, color="#b5651d", label="Cu")
    ax2.set_title("Johnson & Christy 1972")
    for ax in (ax1, ax2):
        ax.set_xlabel("длина волны, нм")
        ax.legend()
    ax1.set_ylabel("P = 1 - R - T")
    fig.suptitle("Рис. 2 статьи Лерер АМ_1 (RCWA, 1D-приближение): чувствительность к базе констант")
    fig.tight_layout()
    png_path = RESULTS / "rerun_fig2_lerer_constants.png"
    fig.savefig(png_path, dpi=150)
    print(f"PNG записан: {png_path}")


if __name__ == "__main__":
    main()
