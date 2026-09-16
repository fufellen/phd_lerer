"""Рисунок по таблице каскада (taper_eme_results.csv): потери перехода против формы острия и длины, в двух мерах -
потери самого перехода (T_fwd, в дБ) и мера трёхмерной модели (T_fwd_eq, относительно перехода без металла сверх
рабочего); плюс потери преобразования (T/T_ад) и ветвь alpha(w) с окном перестройки.

Запуск: python plot_taper_eme.py <taper_eme_results.csv> <taper_branch.csv> <папка вывода>
"""
import csv
import io
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

RU = {"linear": "линейный", "quad": "парабола, острый", "sqrt": "парабола, тупой", "ellipse": "эллипс, тупой",
      "ellipsec": "эллипс, острый", "expo": "экспонента", "gauss": "гаусс", "rcos": "приподн. косинус",
      "klop": "Клопфенштейн", "fastwin": "окно за 5 %", "slowwin": "окно за 50 %", "step": "скачок при z = L",
      "cubic": "куб", "butt": "встык, без участка", "table": "подобранный"}
ORDER = ["butt", "step", "ellipsec", "expo", "cubic", "quad", "gauss", "rcos", "klop", "linear", "fastwin", "ellipse", "sqrt", "slowwin"]


def main():
    eme_path, branch_path, dst = sys.argv[1], sys.argv[2], sys.argv[3]
    rows = list(csv.DictReader(io.open(eme_path, encoding="utf-8")))
    Ls = sorted(set(float(r["L_um"]) for r in rows))
    shapes = [s for s in ORDER if any(r["shape"] == s for r in rows)]
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.2))
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(Ls)))
    x = np.arange(len(shapes))
    w = 0.8 / len(Ls)
    for k, L in enumerate(Ls):
        d = {r["shape"]: r for r in rows if float(r["L_um"]) == L}
        loss = [-10 * math.log10(float(d[s]["T_fwd"])) if s in d else np.nan for s in shapes]
        loss_eq = [-10 * math.log10(float(d[s]["T_fwd_eq"])) if s in d else np.nan for s in shapes]
        conv = [float(d[s]["conv_fwd_dB"]) if s in d else np.nan for s in shapes]
        axes[0].bar(x + (k - (len(Ls) - 1) / 2) * w, loss, w, color=colors[k], label="L = %.2f мкм" % L)
        axes[1].bar(x + (k - (len(Ls) - 1) / 2) * w, loss_eq, w, color=colors[k], label="L = %.2f мкм" % L)
        axes[2].bar(x + (k - (len(Ls) - 1) / 2) * w, conv, w, color=colors[k], label="L = %.2f мкм" % L)
    for ax, title, ylab in ((axes[0], "потери перехода как такового", "−10 lg T, дБ (доля мощности в рабочей моде в конце перехода)"),
                            (axes[1], "мера трёхмерной модели: относительно перехода без металла сверх рабочего", "−10 lg (T·exp(α₄₅₀L)), дБ; ниже нуля — меньше металла, чем в рабочей секции"),
                            (axes[2], "потери преобразования: T относительно exp(−∫α dz)", "−10 lg (T/T_ад), дБ")):
        ax.set_xticks(x)
        ax.set_xticklabels([RU[s] for s in shapes], rotation=35, ha="right", fontsize=8)
        ax.set_title(title, fontsize=9)
        ax.set_ylabel(ylab, fontsize=8)
        ax.grid(alpha=0.3, axis="y")
        ax.axhline(0, color="k", lw=0.6)
        ax.legend(fontsize=7)
    fig.suptitle("Каскад местных мод (femwell, 100 мод, 25 сечений): переход по ширине золота против формы острия и длины", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(dst, "taper_eme_shapes.png"), dpi=150)

    # ветвь alpha(w)
    br = list(csv.DictReader(io.open(branch_path, encoding="utf-8")))
    wv = np.array([float(r["w_au_nm"]) for r in br if int(r["w_au_nm"]) >= 0])
    al = np.array([float(r["alpha_1um"]) for r in br if int(r["w_au_nm"]) >= 0])
    re = np.array([float(r["re_neff"]) for r in br if int(r["w_au_nm"]) >= 0])
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    ax.plot(wv, al, "o-", color="#d62728", label="коэффициент затухания α, мкм⁻¹")
    ax.set_xlabel("ширина золота под гребнем 450 нм, нм")
    ax.set_ylabel("α местной моды, мкм⁻¹", color="#d62728")
    ax.axvspan(210, 290, color="gray", alpha=0.2, label="окно перестройки 210…290 нм")
    ax2 = ax.twinx()
    ax2.plot(wv, re, "s--", color="#1f77b4", label="Re n_eff")
    ax2.set_ylabel("Re n_eff местной моды", color="#1f77b4")
    ax.set_title("Прослеженная ветвь местной моды сечения перехода (femwell, продолжение по ближайшему корню)", fontsize=9)
    ax.grid(alpha=0.3)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(os.path.join(dst, "taper_branch.png"), dpi=150)
    print("taper_eme_shapes.png, taper_branch.png")


if __name__ == "__main__":
    main()
