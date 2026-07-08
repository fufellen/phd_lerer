"""
Валидация ветки A (феноменологическая chi^(3)) против эксперимента
Ryasnyanskiy 2007/2009: предсказанная эффективная восприимчивость композита
chi^(3)_eff = p * |f|^4 * |chi_m| сравнивается с ИЗМЕРЕННОЙ composite-величиной
|chi^(3)_eff| для трёх матриц (Al2O3 / SiO2 / ZnO) на 532 нм, 7 нс.

Зачем этот скрипт. Флагманская фигура P(lambda, I)
(grating_2d_rcwa/rcwa_nonlinear_au_absorption.py) - это ВЫХОД модели для
гипотетической решётки, которую никто не измерял, поэтому нанести на неё
экспериментальные точки нельзя. Сравнение с экспериментом живёт на уровне
chi^(3)_eff композита: именно её измерил Z-scan Ryasnyanskiy. Здесь и строится
"модель vs эксперимент" с литературными точками.

Модель (дилютный предел, ветка A, конвенция Лерера e^{+iwt}):
    f          = 3*eps_h / (eps_p + 2*eps_h)      (квазистатический local-field)
    |chi_eff|  = p * |f|^4 * |chi_m|              (Hache 1988 Eq.(1);
                                                   Ryasnyanskiy Eq.(1)-(2);
                                                   Sipe-Boyd 1992)
eps_p = eps_au(532) из той же таблицы, что и линейный код (composite_ema.py);
|f| и, значит, |f|^4 инвариантны к знаку Im(eps), поэтому конвенция
eps=eps'-i*eps'' здесь роли не играет (сравниваются модули).

ОСТОРОЖНО - полу-круговой характер сверки: Ryasnyanskiy извлекли |chi_m| ИЗ
измеренной |chi_eff| той же local-field теорией. Поэтому это проверка
СОГЛАСОВАННОСТИ нашей реализации f / Максвелла-Гарнетта с их извлечением, на
НЕЗАВИСИМЫХ оптических константах Au (наша таблица != их таблица), а не полностью
независимое измерение chi_eff. Тождеством это не является: разные eps_au(532) и
|f|^4 => расхождение показывает чувствительность множителя усиления |f|^4 к
источнику констант у резонанса. Совпадение по всем трём матрицам в пределах
фактора ~2 (ratio 0.80 / 0.75 / 1.10) подтверждает, что |f|^4 в ветке A
устойчив, а не подогнан.

Показатели преломления матриц - стандартные литературные при 532 нм
(Al2O3 1.77 - совпадает с host композита в статье Лерера и с committed
self_test nonlinear_au_sweep.py; SiO2 1.46; ZnO 2.00). Точные ratio для
SiO2/ZnO несут неопределённость n матрицы; опорная точка - Al2O3 (ratio 0.80).

Числа - из сверенной по PDF таблицы Obsidian-заметки "Нелинейность наночастиц
Au для EMA и расчетов композита", раздел "Сверка benchmarks по PDF
(2026-07-06)", строка Ryasnyanskiy 2007=2009.

Статус: validated (order-of-magnitude / consistency-check).
Выход: composite_ema/results/ryasn_chi3_validation.{csv,png} + печать самопроверки.

Запуск (из корня репозитория):
    python composite_ema/scripts/ryasn_chi3_validation.py
"""

from __future__ import annotations

import csv
from dataclasses import dataclass

from composite_ema import RESULTS_DIR, eps_au
from nonlinear_au_sweep import local_field_factor

LAM_NM = 532.0  # длина волны Z-scan Ryasnyanskiy


@dataclass(frozen=True)
class RyasnPoint:
    matrix: str
    n_host: float           # показатель преломления матрицы при 532 нм (стандартный)
    p_fill: float           # объёмная доля Au ПНЧ
    chi_m_esu: float        # МОДУЛЬ intrinsic chi^(3) частицы (Ryasnyanskiy)
    re_chi_eff_esu: float   # Re chi^(3)_eff композита (измерено)
    im_chi_eff_esu: float   # Im chi^(3)_eff композита (измерено; 0 = не наблюдалось)

    @property
    def meas_abs_chi_eff(self) -> float:
        return abs(complex(self.re_chi_eff_esu, self.im_chi_eff_esu))


# Сверенная по PDF таблица (Obsidian: "Сверка benchmarks по PDF (2026-07-06)").
POINTS = [
    RyasnPoint("Al2O3", 1.77, 0.0800,  6.25e-8,  5.03e-7, -5.33e-7),
    RyasnPoint("SiO2",  1.46, 0.0800,  2.70e-8,  1.38e-7, -0.067e-7),
    RyasnPoint("ZnO",   2.00, 0.0868,  2.46e-8, -1.47e-7,  0.0),
]


def predict_abs_chi_eff(pt: RyasnPoint) -> tuple[float, float]:
    """Возвращает (|chi_eff|_pred в esu, |f|^4) для точки."""
    eps_h = complex(pt.n_host ** 2, 0.0)
    f = local_field_factor(eps_au(LAM_NM), eps_h)
    f4 = abs(f) ** 4
    return pt.p_fill * f4 * pt.chi_m_esu, f4


def self_test() -> None:
    # Опорная точка Al2O3 должна совпасть с committed self_test
    # nonlinear_au_sweep.py (ratio 0.80).
    al = POINTS[0]
    pred, _ = predict_abs_chi_eff(al)
    ratio = pred / al.meas_abs_chi_eff
    assert abs(ratio - 0.80) < 0.03, ratio
    print(f"[self_test] OK: Al2O3 ratio pred/meas = {ratio:.2f} "
          f"(совпадает с self_test nonlinear_au_sweep.py = 0.80).")
    # Все три матрицы - в пределах фактора 2 от диагонали y=x.
    for pt in POINTS:
        pred, _ = predict_abs_chi_eff(pt)
        r = pred / pt.meas_abs_chi_eff
        assert 0.5 < r < 2.0, (pt.matrix, r)
    print("[self_test] OK: все три матрицы в пределах фактора 2 от y=x.")


def run() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    rows = []
    for pt in POINTS:
        pred, f4 = predict_abs_chi_eff(pt)
        meas = pt.meas_abs_chi_eff
        rows.append({
            "matrix": pt.matrix,
            "n_host": pt.n_host,
            "p_fill": pt.p_fill,
            "chi_m_esu": pt.chi_m_esu,
            "abs_f4": f4,
            "pred_abs_chi_eff_esu": pred,
            "meas_abs_chi_eff_esu": meas,
            "ratio_pred_over_meas": pred / meas,
        })

    csv_path = RESULTS_DIR / "ryasn_chi3_validation.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"[run] CSV записан: {csv_path}")
    for r in rows:
        print(f"[run]   {r['matrix']}: pred={r['pred_abs_chi_eff_esu']:.3e} esu  "
              f"meas={r['meas_abs_chi_eff_esu']:.3e} esu  "
              f"ratio={r['ratio_pred_over_meas']:.2f}  (|f|^4={r['abs_f4']:.1f})")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    fig, ax = plt.subplots(figsize=(6.6, 6.2))
    lo, hi = 5e-8, 2e-6
    xx = np.array([lo, hi])
    ax.plot(xx, xx, color="gray", lw=1.2, label="идеальное совпадение (y=x)")
    ax.fill_between(xx, xx / 2.0, xx * 2.0, color="gray", alpha=0.12, label="фактор 2")

    colors = {"Al2O3": "#c02b2b", "SiO2": "#2f6f8f", "ZnO": "#2a8f4f"}
    for r in rows:
        x = r["meas_abs_chi_eff_esu"]
        y = r["pred_abs_chi_eff_esu"]
        ax.scatter([x], [y], s=95, color=colors[r["matrix"]], zorder=5,
                   edgecolor="black", linewidth=0.6)
        ax.annotate(f"Au:{r['matrix']}\n×{r['ratio_pred_over_meas']:.2f}",
                    (x, y), textcoords="offset points", xytext=(10, -6), fontsize=9)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal")
    ax.set_xlabel(r"измеренная $|\chi^{(3)}_\mathrm{eff}|$, esu  (Ryasnyanskiy 2007)")
    ax.set_ylabel(r"модель  $p\,|f|^4\,|\chi^{(3)}_m|$, esu")
    ax.set_title("Валидация ветки A: composite $\\chi^{(3)}$ модель vs эксперимент\n"
                 "Ryasnyanskiy 2007 (Au ПНЧ, 532 нм, 7 нс, 3 матрицы)", fontsize=10)
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(True, which="both", alpha=0.25)
    fig.tight_layout()

    png_path = RESULTS_DIR / "ryasn_chi3_validation.png"
    fig.savefig(png_path, dpi=150)
    print(f"[run] PNG записан: {png_path}")


if __name__ == "__main__":
    self_test()
    run()
