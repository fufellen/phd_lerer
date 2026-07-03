"""
Стадия 4 (продолжение): геометрия рис. 2 статьи [[Лерер АМ_1]] - теперь с
НАСТОЯЩЕЙ 2D-решёткой цилиндров (rcwa_2d.py, crossed grating), а не
1D-полосками (reproduce_fig2_au_vs_cu.py).

Та же геометрия из подписи к рис. 2 статьи: периоды 500 нм (квадратная
решётка), диаметр цилиндров 400 нм, высота 100 нм; под ними сплошной слой
с ПНЧ h=100 нм; подложка n=1.45; матрица n=1.77, C=10%.

Ключевое отличие от 1D: заполняющая доля площади у цилиндров
pi*(D/2)^2/Λ^2 ~ 0.503, а не D/Λ = 0.8 как у полосок - в 2D в ячейке
меньше поглощающего материала. Плюс появляются косые дифракционные
порядки (m,n) и нет искусственной анизотропии полосок.

Цель: проверить, воспроизводится ли УФ-сравнение Au vs Cu из статьи
(там Au > Cu в УФ) в истинной 2D-геометрии - в 1D-приближении оно
НЕ воспроизводилось (Cu > Au). Это была одна из трёх гипотез о причине
расхождения в отчёте о воспроизведении.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "composite_ema" / "scripts"))

from rcwa_2d import fourier_coeffs_cylinder, convolution_matrix_2d, solve_rcwa_2d  # noqa: E402
import composite_ema as ema  # noqa: E402

N_HARM = 5  # (2N+1)^2 = 121 гармоника; P сходится уже при N=4 (см. validate_rcwa_2d)


def run(material_n_mat: int, material_name: str, N: int = N_HARM):
    n_host = 1.77
    C = 0.10
    n_sub = 1.45
    period = 500.0
    diameter = 400.0
    h_cyl = 100.0
    h_layer = 100.0
    eps_in = 1.0 + 1e-9j  # воздух сверху; крошечная мнимая часть - от аномалий Вуда (см. 1D-скрипт)
    eps_out = n_sub ** 2 + 1e-9j

    lo, hi = ema.MATERIAL_RANGE_NM[material_n_mat]
    lam_grid = np.arange(max(210.0, lo), min(900.0, hi), 5.0)

    R_list, T_list, P_list = [], [], []
    t0 = time.perf_counter()
    for i, lam in enumerate(lam_grid):
        eps_h = complex(n_host ** 2, 0.0)
        eps_p = ema.libmat(material_n_mat, float(lam))
        eps_comp = ema.composite_from_eps(eps_h, eps_p, C)   # конвенция Лерера eps'-i*eps''
        eps_comp_std = eps_comp.conjugate()                   # стандартная фотонная конвенция

        a = fourier_coeffs_cylinder(eps_comp_std, 1.0 + 0j, diameter / period, N)
        E_cyl = convolution_matrix_2d(a, N)
        layers = [(E_cyl, h_cyl), (eps_comp_std, h_layer)]   # скаляр -> однородный слой (аналитический базис)

        R_mn, T_mn = solve_rcwa_2d(layers, eps_in, eps_out, N, period, float(lam), pol="y")
        R = float(R_mn.sum().real)
        Tt = float(T_mn.sum().real)
        R_list.append(R)
        T_list.append(Tt)
        P_list.append(1.0 - R - Tt)
        if i % 20 == 0:
            print(f"  [{material_name}] {lam:.0f} нм ({i+1}/{len(lam_grid)}), "
                  f"P={P_list[-1]:.3f}, {time.perf_counter()-t0:.0f} с")

    R_arr, T_arr, P_arr = np.array(R_list), np.array(T_list), np.array(P_list)
    print(f"[{material_name}] диапазон {lam_grid[0]:.0f}-{lam_grid[-1]:.0f} нм, "
          f"максимум P={P_arr.max():.3f} на {lam_grid[np.argmax(P_arr)]:.0f} нм "
          f"({time.perf_counter()-t0:.0f} с)")
    return lam_grid, R_arr, T_arr, P_arr


def main() -> None:
    print(f"2D RCWA (цилиндры), N={N_HARM} -> {(2*N_HARM+1)**2} гармоник")
    lam_au, R_au, T_au, P_au = run(4, "Au")
    lam_cu, R_cu, T_cu, P_cu = run(1, "Cu")

    assert np.all(P_au > -1e-6) and np.all(P_cu > -1e-6), "P<0 - нефизично, проверить конвенцию знака"
    print("[OK] P>=0 всюду (пассивные материалы).")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.5, 5))
    ax.plot(lam_au, P_au, label="Au (2D цилиндры)", color="#d4a017")
    ax.plot(lam_cu, P_cu, label="Cu (2D цилиндры)", color="#b5651d")
    ax.set_xlabel("wavelength, nm")
    ax.set_ylabel("P (доля поглощённой мощности)")
    ax.set_title("RCWA 2D (цилиндры) - рис. 2 статьи: период 500нм, диаметр 400нм, C=10%")
    ax.legend()
    fig.tight_layout()
    out_dir = Path(__file__).resolve().parents[1] / "results"
    out_dir.mkdir(parents=True, exist_ok=True)
    png_path = out_dir / "reproduce_fig2_au_vs_cu_2d.png"
    fig.savefig(png_path, dpi=150)
    print(f"PNG записан: {png_path}")

    csv_path = out_dir / "reproduce_fig2_au_vs_cu_2d.csv"
    with csv_path.open("w", encoding="utf-8") as fh:
        fh.write("lambda_nm,R_Au,T_Au,P_Au,lambda_nm_Cu,R_Cu,T_Cu,P_Cu\n")
        n = max(len(lam_au), len(lam_cu))
        for i in range(n):
            row = [
                str(lam_au[i]) if i < len(lam_au) else "",
                str(R_au[i]) if i < len(lam_au) else "",
                str(T_au[i]) if i < len(lam_au) else "",
                str(P_au[i]) if i < len(lam_au) else "",
                str(lam_cu[i]) if i < len(lam_cu) else "",
                str(R_cu[i]) if i < len(lam_cu) else "",
                str(T_cu[i]) if i < len(lam_cu) else "",
                str(P_cu[i]) if i < len(lam_cu) else "",
            ]
            fh.write(",".join(row) + "\n")
    print(f"CSV записан: {csv_path}")

    print("\n--- Оценка (2D цилиндры) ---")
    peak_au_lam = lam_au[np.argmax(P_au)]
    peak_cu_lam = lam_cu[np.argmax(P_cu)]
    print(f"Резонанс Au (P max): {peak_au_lam:.0f} нм; Cu (P max): {peak_cu_lam:.0f} нм")
    print(f"Условие Фрёлиха (composite_ema): Au ~540 нм, Cu ~528 нм")
    print(f"P_max Au = {P_au.max():.3f}, P_max Cu = {P_cu.max():.3f}")

    uv_mask_au = lam_au <= 400.0
    uv_mask_cu = lam_cu <= 400.0
    if uv_mask_au.any() and uv_mask_cu.any():
        p_au_uv = P_au[uv_mask_au].mean()
        p_cu_uv = P_cu[uv_mask_cu].mean()
        print(f"УФ (<=400 нм), среднее P: Au = {p_au_uv:.3f}, Cu = {p_cu_uv:.3f}")
        if p_au_uv > p_cu_uv:
            print("2D-геометрия СОГЛАСУЕТСЯ с рис. 2 статьи: Au > Cu в УФ "
                  "(в 1D-приближении было наоборот - гипотеза о роли геометрии подтверждается).")
        else:
            print("2D-геометрия НЕ снимает расхождение: Cu > Au в УФ и с цилиндрами. "
                  "Остаётся гипотеза о разных первоисточниках оптических констант "
                  "(luxpop у Лерера против Johnson&Christy здесь).")


if __name__ == "__main__":
    main()
