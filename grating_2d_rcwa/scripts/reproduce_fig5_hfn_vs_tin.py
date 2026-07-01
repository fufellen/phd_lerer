"""
Стадия 3, рис. 5 статьи [[Лерер АМ_1]]: HfN vs TiN, та же геометрия и то же
1D-приближение решётки, что в reproduce_fig2_au_vs_cu.py (см. его докстринг
и README.md "Упрощение геометрии" - применимо один в один, только материал
меняется).

Цель проверки: статья утверждает (см. "Заключение" и подпись рис. 5), что
TiN даёт более широкополосное поглощение, чем HfN. Тот же качественный
вывод уже был получен на уровне composite_ema (formula (1) без решётки,
см. composite_ema/README.md - полоса TiN 346нм vs HfN 208нм). Здесь
проверяется, сохраняется ли это соотношение после добавления настоящей
дифракции на решётке (не только эффективная eps слоя).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "composite_ema" / "scripts"))

from rcwa_1d_te import fourier_coeffs_strip, toeplitz_from_coeffs, solve_rcwa_te  # noqa: E402
import composite_ema as ema  # noqa: E402


def run(material_n_mat: int, material_name: str, N: int = 8) -> tuple[np.ndarray, np.ndarray]:
    n_host = 1.77
    C = 0.10
    n_sub = 1.45
    period = 500.0
    diameter = 400.0
    fill = diameter / period
    h_cyl = 100.0
    h_layer = 100.0
    eps_in = 1.0 + 1e-9j
    eps_out = n_sub ** 2 + 1e-9j

    lo, hi = ema.MATERIAL_RANGE_NM[material_n_mat]
    lam_grid = np.arange(max(350.0, lo), min(900.0, hi), 5.0)

    P_list = []
    for lam in lam_grid:
        eps_h = complex(n_host ** 2, 0.0)
        eps_p = ema.libmat(material_n_mat, float(lam))
        eps_comp = ema.composite_from_eps(eps_h, eps_p, C).conjugate()
        coeffs = fourier_coeffs_strip(eps_comp, 1.0 + 0j, fill, N)
        T_cyl = toeplitz_from_coeffs(coeffs)
        coeffs_flat = fourier_coeffs_strip(eps_comp, eps_comp, 0.0, N)
        T_flat = toeplitz_from_coeffs(coeffs_flat)
        layers = [(T_cyl, h_cyl), (T_flat, h_layer)]
        R_harm, T_harm = solve_rcwa_te(layers, eps_in, eps_out, N, period, float(lam), theta_i_deg=0.0)
        R = float(np.sum(R_harm).real)
        Tt = float(np.sum(T_harm).real)
        P_list.append(1.0 - R - Tt)

    P_arr = np.array(P_list)
    half = P_arr.max() / 2.0
    above = lam_grid[P_arr >= half]
    bw = (above.max() - above.min()) if len(above) else 0.0
    print(f"[{material_name}] P_max={P_arr.max():.3f} на {lam_grid[np.argmax(P_arr)]:.0f} нм, "
          f"полоса half-max={bw:.0f} нм ({above.min():.0f}-{above.max():.0f})")
    return lam_grid, P_arr


def main() -> None:
    lam_hfn, P_hfn = run(5, "HfN")
    lam_tin, P_tin = run(7, "TiN")

    half_hfn = P_hfn.max() / 2.0
    above_hfn = lam_hfn[P_hfn >= half_hfn]
    bw_hfn = above_hfn.max() - above_hfn.min()
    half_tin = P_tin.max() / 2.0
    above_tin = lam_tin[P_tin >= half_tin]
    bw_tin = above_tin.max() - above_tin.min()

    print(f"\nПолоса HfN: {bw_hfn:.0f} нм, полоса TiN: {bw_tin:.0f} нм")
    if bw_tin > bw_hfn:
        print("Качественно согласуется со статьёй: полоса TiN шире полосы HfN (подтверждено и на уровне полной решётки, не только composite_ema).")
    else:
        print("РАСХОЖДЕНИЕ со статьёй: ожидалась более широкая полоса у TiN.")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.5, 5))
    ax.plot(lam_hfn, P_hfn, label="HfN, P=1-R-T", color="#8b3a62")
    ax.plot(lam_tin, P_tin, label="TiN, P=1-R-T", color="#2f6f8f")
    ax.set_xlabel("wavelength, nm")
    ax.set_ylabel("P (доля поглощённой мощности)")
    ax.set_title("RCWA (1D-приближение решётки) - рис. 5 статьи: HfN vs TiN, период 500нм, диаметр 400нм")
    ax.legend()
    fig.tight_layout()
    out_dir = Path(__file__).resolve().parents[1] / "results"
    out_dir.mkdir(parents=True, exist_ok=True)
    png_path = out_dir / "reproduce_fig5_hfn_vs_tin.png"
    fig.savefig(png_path, dpi=150)
    print(f"PNG записан: {png_path}")

    csv_path = out_dir / "reproduce_fig5_hfn_vs_tin.csv"
    with csv_path.open("w", encoding="utf-8") as fh:
        fh.write("lambda_nm_HfN,P_HfN,lambda_nm_TiN,P_TiN\n")
        n = max(len(lam_hfn), len(lam_tin))
        for i in range(n):
            row = [
                str(lam_hfn[i]) if i < len(lam_hfn) else "",
                str(P_hfn[i]) if i < len(lam_hfn) else "",
                str(lam_tin[i]) if i < len(lam_tin) else "",
                str(P_tin[i]) if i < len(lam_tin) else "",
            ]
            fh.write(",".join(row) + "\n")
    print(f"CSV записан: {csv_path}")


if __name__ == "__main__":
    main()
