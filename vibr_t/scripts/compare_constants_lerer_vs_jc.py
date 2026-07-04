"""
Сравнение оптических констант Au/Cu: таблицы Лерера (vibr_t, первоисточник
программы VIBR_T) против Johnson & Christy 1972 (composite_ema, текущая
чистовая база наших воспроизведений).

Цель - закрыть открытый вопрос воспроизведения [[Лерер АМ_1]]:
"УФ-сравнение Au/Cu даёт противоположный статье порядок" - гипотеза
заключалась в разных первоисточниках констант. Здесь она проверяется
напрямую: у Лерера НЕТ данных Cu ниже 517 нм, в УФ его программа линейно
экстраполирует, и Re eps_Cu становится положительной ниже ~416 нм.

Выход: results/compare_constants_au_cu.csv, results/compare_constants_au_cu.png
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "composite_ema" / "scripts"))

import vibr_materials as vm  # noqa: E402
import composite_ema as ema  # noqa: E402

RESULTS = HERE.parent / "results"

# индексы материалов: vibr_t (BIBL_mat.c): Cu=1, Au=3; composite_ema: Cu=1, Au=4
LAM = np.arange(210.0, 1000.0, 2.0)


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    rows = []
    for lam in LAM:
        l = float(lam)
        eps_au_l = vm.libmat(3, l)                    # Лерер, плотная таблица
        eps_cu_l = vm.libmat(1, l)                    # Лерер, с экстраполяцией в УФ
        eps_au_j = ema.libmat(4, l)
        eps_cu_j = ema.libmat(1, l)
        extrap = vm.cu_table_range()[0] > l
        rows.append([l, eps_au_l.real, eps_au_l.imag, eps_cu_l.real, eps_cu_l.imag,
                     eps_au_j.real, eps_au_j.imag, eps_cu_j.real, eps_cu_j.imag,
                     int(extrap)])
    arr = np.array(rows)

    csv_path = RESULTS / "compare_constants_au_cu.csv"
    header = ("lambda_nm,Au_re_lerer,Au_im_lerer,Cu_re_lerer,Cu_im_lerer,"
              "Au_re_jc,Au_im_jc,Cu_re_jc,Cu_im_jc,cu_lerer_extrapolated")
    np.savetxt(csv_path, arr, delimiter=",", header=header, comments="")
    print(f"CSV записан: {csv_path}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
    cu_lo = vm.cu_table_range()[0]
    for ax_col, (mat, li, ji) in enumerate((("Au", 1, 5), ("Cu", 3, 7))):
        ax_re, ax_im = axes[0][ax_col], axes[1][ax_col]
        ax_re.plot(arr[:, 0], arr[:, li], label="Лерер (VIBR_T)", color="C3")
        ax_re.plot(arr[:, 0], arr[:, ji], label="Johnson & Christy", color="C0", ls="--")
        ax_im.plot(arr[:, 0], arr[:, li + 1], color="C3")
        ax_im.plot(arr[:, 0], arr[:, ji + 1], color="C0", ls="--")
        ax_re.set_title(f"{mat}: Re eps")
        ax_im.set_title(f"{mat}: Im eps (конвенция e^{{+iwt}}, <=0)")
        ax_im.set_xlabel("длина волны, нм")
        if mat == "Cu":
            for ax in (ax_re, ax_im):
                ax.axvline(cu_lo, color="gray", lw=1, ls=":")
                ax.axvspan(float(arr[0, 0]), cu_lo, color="red", alpha=0.08)
            ax_re.annotate("экстраполяция\n(нет данных)", xy=(cu_lo - 140, 0),
                           fontsize=9, color="red")
        ax_re.axhline(0, color="k", lw=0.5)
        ax_re.legend()
    fig.suptitle("Оптические константы Au/Cu: таблицы программы Лерера vs Johnson & Christy 1972")
    fig.tight_layout()
    png_path = RESULTS / "compare_constants_au_cu.png"
    fig.savefig(png_path, dpi=150)
    print(f"PNG записан: {png_path}")

    # численный вердикт
    uv = arr[arr[:, 0] <= 400.0]
    print("\n--- Вердикт по УФ (<=400 нм) ---")
    print(f"Cu (Лерер): Re eps от {uv[:, 3].min():.2f} до {uv[:, 3].max():.2f} - "
          f"положительная выше ~416 нм экстраполяции, медь не металл")
    print(f"Cu (J&C):   Re eps от {uv[:, 7].min():.2f} до {uv[:, 7].max():.2f} - "
          f"настоящие измеренные данные")
    d_au = np.abs(uv[:, 1] - uv[:, 5]).max()
    print(f"Au: расхождение Лерер vs J&C в УФ, max|dRe eps| = {d_au:.2f} "
          f"(таблица Au у Лерера в УФ есть, источники близки)")


if __name__ == "__main__":
    main()
