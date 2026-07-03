"""
Стадия 3 (первая попытка): геометрия рис. 2 статьи [[Лерер АМ_1]] - через
RCWA (1D-приближение решётки, см. rcwa_1d_te.py) + материалы composite_ema.

Параметры из подписи к рис. 2 статьи: "Периоды 500 нм, диаметр цилиндров
400 нм" + из текста "Результаты": подложка n=1.45, диэлектрический слой с
ПНЧ h=100нм n=1.77 C=10%, цилиндры h=100нм n=1.77 C=10%.

Упрощение: реальные цилиндры (2D-периодические, круглые) заменены на 1D
полоски той же заполняющей доли (fill = diameter/period = 0.8) - см.
README.md про то, почему 1D, а не 2D. Оптическая ось: TE (нормальное
падение, Ey поперёк полосок).

Не ожидается количественного совпадения с рис. 2 (другая геометрия
рассеяния, 1D vs 2D) - цель: качественное поведение (резонанс, Au vs Cu),
см. README.md "Цель задачи".
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "composite_ema" / "scripts"))

from rcwa_1d_te import fourier_coeffs_strip, toeplitz_from_coeffs, solve_rcwa_te  # noqa: E402
import composite_ema as ema  # noqa: E402


def eps_layer_and_grating(n_mat: int, lam_nm: float, n_host: float, C: float) -> complex:
    """eps_eff композита (host+ПНЧ) для материала n_mat при данной lambda -
    используется и для сплошного слоя, и для материала внутри полосок."""
    eps_h = complex(n_host ** 2, 0.0)
    eps_p = ema.libmat(n_mat, lam_nm)
    return ema.composite_from_eps(eps_h, eps_p, C)


def run(material_n_mat: int, material_name: str, N: int = 6) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n_host = 1.77
    C = 0.10
    n_sub = 1.45
    period = 500.0
    diameter = 400.0
    fill = diameter / period
    h_cyl = 100.0
    h_layer = 100.0
    # Крошечная мнимая добавка нужна только для того, чтобы численно избежать
    # точного попадания на аномалию Вуда (kz=0 на какой-то гармонике,
    # np.linalg.solve падает на точно вырожденной матрице) при сканировании
    # по сетке длин волн - физику не меняет (беспотерьный воздух/подложка).
    eps_in = 1.0 + 1e-9j  # воздух сверху
    eps_out = n_sub ** 2 + 1e-9j

    lo, hi = ema.MATERIAL_RANGE_NM[material_n_mat]
    lam_grid = np.arange(max(210.0, lo), min(900.0, hi), 5.0)

    R_list, T_list, P_list = [], [], []
    for lam in lam_grid:
        eps_comp = eps_layer_and_grating(material_n_mat, float(lam), n_host, C)  # composite_ema convention (eps'-i*eps'')
        eps_comp_std = eps_comp.conjugate()  # -> стандартная фотонная конвенция для tmm/rcwa (см. grating_2d_rcwa/README.md)

        # слой цилиндров: полоски из eps_comp_std в фоне воздуха (eps=1)
        coeffs = fourier_coeffs_strip(eps_comp_std, 1.0 + 0j, fill, N)
        T_cyl = toeplitz_from_coeffs(coeffs)

        # сплошной диэлектрический слой с ПНЧ (без решётки - только DC член)
        coeffs_flat = fourier_coeffs_strip(eps_comp_std, eps_comp_std, 0.0, N)
        T_flat = toeplitz_from_coeffs(coeffs_flat)

        layers = [(T_cyl, h_cyl), (T_flat, h_layer)]
        R_harm, T_harm = solve_rcwa_te(layers, eps_in, eps_out, N, period, float(lam), theta_i_deg=0.0)
        R = float(np.sum(R_harm).real)
        Tt = float(np.sum(T_harm).real)
        R_list.append(R)
        T_list.append(Tt)
        P_list.append(1.0 - R - Tt)

    R_arr, T_arr, P_arr = np.array(R_list), np.array(T_list), np.array(P_list)
    print(f"[{material_name}] диапазон {lam_grid[0]:.0f}-{lam_grid[-1]:.0f} нм, "
          f"максимум P={P_arr.max():.3f} на {lam_grid[np.argmax(P_arr)]:.0f} нм")
    return lam_grid, R_arr, T_arr, P_arr


def main() -> None:
    lam_au, R_au, T_au, P_au = run(4, "Au")
    lam_cu, R_cu, T_cu, P_cu = run(1, "Cu")

    max_imbalance = max(np.max(np.abs(R_au + T_au + P_au - 1.0)), np.max(np.abs(R_cu + T_cu + P_cu - 1.0)))
    print(f"[energy] максимум |R+T+P-1| (должно быть ~0, здесь среды с потерями, R+T+P=1 всегда по построению): {max_imbalance:.2e}")
    assert np.all(P_au > -1e-6) and np.all(P_cu > -1e-6), "P<0 - нефизичное усиление, проверить конвенцию знака"
    print("[OK] P>=0 всюду (пассивные материалы, поглощение физично).")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.5, 5))
    ax.plot(lam_au, P_au, label="Au, P=1-R-T", color="#d4a017")
    ax.plot(lam_cu, P_cu, label="Cu, P=1-R-T", color="#b5651d")
    ax.set_xlabel("wavelength, nm")
    ax.set_ylabel("P (доля поглощённой мощности)")
    ax.set_title("RCWA (1D-приближение решётки) - рис. 2 статьи: период 500нм, диаметр 400нм, C=10%")
    ax.legend()
    fig.tight_layout()
    out_dir = Path(__file__).resolve().parents[1] / "results"
    out_dir.mkdir(parents=True, exist_ok=True)
    png_path = out_dir / "reproduce_fig2_au_vs_cu.png"
    fig.savefig(png_path, dpi=150)
    print(f"PNG записан: {png_path}")

    csv_path = out_dir / "reproduce_fig2_au_vs_cu.csv"
    with csv_path.open("w", encoding="utf-8") as fh:
        fh.write("lambda_nm,R_Au,T_Au,P_Au,lambda_nm_Cu,R_Cu,T_Cu,P_Cu\n")
        n = max(len(lam_au), len(lam_cu))
        for i in range(n):
            row = []
            row.append(str(lam_au[i]) if i < len(lam_au) else "")
            row.append(str(R_au[i]) if i < len(lam_au) else "")
            row.append(str(T_au[i]) if i < len(lam_au) else "")
            row.append(str(P_au[i]) if i < len(lam_au) else "")
            row.append(str(lam_cu[i]) if i < len(lam_cu) else "")
            row.append(str(R_cu[i]) if i < len(lam_cu) else "")
            row.append(str(T_cu[i]) if i < len(lam_cu) else "")
            row.append(str(P_cu[i]) if i < len(lam_cu) else "")
            fh.write(",".join(row) + "\n")
    print(f"CSV записан: {csv_path}")

    # Оценка по цели задачи (см. README.md "Цель задачи")
    print("\n--- Оценка по цели задачи ---")
    peak_au_lam = lam_au[np.argmax(P_au)]
    peak_cu_lam = lam_cu[np.argmax(P_cu)]
    print(f"Резонанс Au (P max): {peak_au_lam:.0f} нм; резонанс Cu (P max): {peak_cu_lam:.0f} нм")
    print(f"Условие Фрёлиха (из composite_ema, Johnson & Christy 1972): Au ~555 нм, Cu ~567 нм")
    print(f"P_max Au = {P_au.max():.3f}, P_max Cu = {P_cu.max():.3f} - цель P>0.8: "
          f"{'ДОСТИГНУТА' if max(P_au.max(), P_cu.max()) > 0.8 else 'НЕ достигнута'}")

    uv_mask_au = lam_au <= 400.0
    uv_mask_cu = lam_cu <= 400.0
    if uv_mask_au.any() and uv_mask_cu.any():
        p_au_uv = P_au[uv_mask_au].mean()
        p_cu_uv = P_cu[uv_mask_cu].mean()
        print(f"УФ (<=400 нм), среднее P: Au = {p_au_uv:.3f}, Cu = {p_cu_uv:.3f}")
        if p_au_uv > p_cu_uv:
            print("Согласуется с рис. 2 статьи: Au поглощает в УФ сильнее, чем Cu.")
        else:
            print("РАСХОЖДЕНИЕ с рис. 2 статьи: в этой модели Cu поглощает в УФ сильнее Au "
                  "(см. обсуждение в итоговой заметке о воспроизведении - возможно, разные "
                  "первоисточники оптических констант Au/Cu, или эффект не воспроизводится "
                  "в 1D-приближении решётки).")


if __name__ == "__main__":
    main()
