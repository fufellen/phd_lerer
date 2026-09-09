"""Рисунок к nonlinear_selfconsistent.py: расцепленное и самосогласованное поглощение нелинейной
метаповерхности и поле в композите. Выход: results/nonlinear_selfconsistent.png."""
from __future__ import annotations

import csv
import io
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import lerer_grcwa as L  # noqa: E402


def main():
    rows = list(csv.DictReader(io.open(L.RESULTS / "nonlinear_selfconsistent.csv", encoding="utf-8")))
    I_all = sorted({float(r["I_wcm2"]) for r in rows})
    lam = np.array(sorted({float(r["lambda_nm"]) for r in rows}))
    get = lambda I, key: np.array([float(r[key]) for r in rows if float(r["I_wcm2"]) == I])
    P_lin = get(I_all[0], "P_lin")
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(I_all)))
    fig, axs = plt.subplots(1, 3, figsize=(15.5, 6.0))
    ax = axs[0]
    ax.plot(lam, P_lin, color="k", lw=2, label="линейный расчёт")
    I = I_all[-1]
    ax.plot(lam, get(I, "P_dec"), color="#d62728", lw=1.8, ls="--", label="расцепленное приближение, I = %.0e Вт/см²" % I)
    ax.plot(lam, get(I, "P_sc"), color="#1f77b4", lw=1.8, label="самосогласованный расчёт, I = %.0e Вт/см²" % I)
    ax.set_ylabel("поглощение P = 1 − R − T")
    ax.set_title("Спектр поглощения решётки цилиндров", fontsize=10)
    ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=1)
    ax = axs[1]
    for I, c in zip(I_all, colors):
        ax.plot(lam, get(I, "dP_dec"), color=c, lw=1.6, ls="--")
        ax.plot(lam, get(I, "dP_sc"), color=c, lw=1.8, label="I = %.0e Вт/см²" % I)
    ax.plot([], [], color="gray", ls="--", label="штрих - расцепленное")
    ax.plot([], [], color="gray", ls="-", label="сплошная - самосогласованное")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_ylabel("изменение поглощения ΔP = P(I) − P(0)")
    ax.set_title("Нелинейный сдвиг поглощения", fontsize=10)
    ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2)
    ax = axs[2]
    ax.plot(lam, get(I_all[0], "ratio_cyl_lin"), color="#ff7f0e", lw=1.8, label="цилиндры, линейное поле")
    ax.plot(lam, get(I_all[0], "ratio_film_lin"), color="#2ca02c", lw=1.8, label="слой, линейное поле")
    ax.plot(lam, get(I_all[-1], "ratio_cyl_sc"), color="#ff7f0e", lw=1.2, ls=":", label="цилиндры, I = %.0e" % I_all[-1])
    ax.plot(lam, get(I_all[-1], "ratio_film_sc"), color="#2ca02c", lw=1.2, ls=":", label="слой, I = %.0e" % I_all[-1])
    ax.axhline(1.0, color="k", lw=1.0, ls="--", label="расцепленное приближение: 1")
    ax.set_ylabel("n_h ⟨|E|²⟩ / |E₀|² в композите")
    ax.set_title("Среднее поле в композите против плоской волны в матрице", fontsize=10)
    ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2)
    for ax in axs:
        ax.set_xlabel("длина волны в свободном пространстве, нм")
        ax.grid(alpha=0.3)
    fig.suptitle("Нелинейная метаповерхность (период 500 нм, цилиндры D = 400 нм, Au 10 % Джонсон-Кристи, χ⁽³⁾ Hache 1988): "
                 "среднеполевой самосогласованный расчёт grcwa", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(L.RESULTS / "nonlinear_selfconsistent.png", dpi=170)
    plt.close(fig)
    print("рисунок записан")


if __name__ == "__main__":
    main()
