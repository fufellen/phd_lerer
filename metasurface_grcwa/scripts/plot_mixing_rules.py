"""Рисунки к mixing_rules_on_grating.py: спектры поглощения решётки цилиндров по правилам смешивания
и схема ячейки. Выход: results/mixing_rules_on_grating.png, results/mixing_rules_scheme.png."""
from __future__ import annotations

import csv
import io
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Rectangle, Ellipse  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import lerer_grcwa as L  # noqa: E402

COL = {"MG": "#1f77b4", "Bruggeman": "#d62728", "MLWA15": "#2ca02c", "Lerer": "#ff7f0e",
       "dilute": "#9467bd", "LL": "#8c564b"}
LAB = {"MG": "Максвелл-Гарнетт (КМ)", "Bruggeman": "Бруггеман", "MLWA15": "МГ + MLWA, R = 15 нм", "Lerer": "формула Лерера",
       "dilute": "первый порядок по C", "LL": "Лоренц-Лоренц"}


def read(path):
    rows = list(csv.DictReader(io.open(path, encoding="utf-8")))
    out = {}
    for r in rows:
        out.setdefault((r["table"], r["formula"]), []).append((float(r["lambda_nm"]), float(r["A"])))
    return {k: np.array(v) for k, v in out.items()}


def main():
    main_rows = read(L.RESULTS / "mixing_rules_on_grating.csv")
    extra = read(L.RESULTS / "mixing_rules_on_grating_extra.csv") if (L.RESULTS / "mixing_rules_on_grating_extra.csv").exists() else {}
    fig, axs = plt.subplots(2, 2, figsize=(12.5, 9.0))
    for ax, table, title in ((axs[0, 0], "lerer", "золото по таблице Лерера (Au.txt)"),
                             (axs[0, 1], "jc", "золото по Джонсону-Кристи")):
        for f in ("MG", "Bruggeman", "MLWA15", "Lerer"):
            d = main_rows[(table, f)]
            ax.plot(d[:, 0], d[:, 1], color=COL[f], lw=2, label=LAB[f])
        for f in ("dilute", "LL"):
            if (table, f) in extra:
                d = extra[(table, f)]
                ax.plot(d[:, 0], d[:, 1], color=COL[f], lw=1.2, ls=":", label=LAB[f])
        ax.set_title(title, fontsize=10.5)
        ax.set_xlabel("длина волны в свободном пространстве, нм")
        ax.set_ylabel("поглощение P = 1 − R − T")
        ax.set_ylim(0, 1.02)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8.5, loc="lower left")
    ax = axs[1, 0]
    for f in ("Bruggeman", "MLWA15", "Lerer"):
        for table, ls in (("lerer", "-"), ("jc", "--")):
            d = main_rows[(table, f)]
            base = main_rows[(table, "MG")]
            ax.plot(d[:, 0], d[:, 1] - base[:, 1], color=COL[f], lw=1.8, ls=ls,
                    label="%s, %s" % (LAB[f], "Au.txt" if table == "lerer" else "Джонсон-Кристи"))
    ax.axhline(0, color="k", lw=0.8)
    ax.set_title("Отличие от Максвелла-Гарнетта (сплошные - Au.txt, штриховые - Джонсон-Кристи)", fontsize=9.5)
    ax.set_xlabel("длина волны в свободном пространстве, нм")
    ax.set_ylabel("P(формула) − P(Максвелл-Гарнетт)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="lower left")
    ax = axs[1, 1]
    labels = {"Au10": "Au 10 %", "Ag10": "Ag 10 %", "Cu10": "Cu 10 %", "Au5+Ag5": "Au 5 % + Ag 5 %",
              "Au5+Cu5": "Au 5 % + Cu 5 %", "COMPOSITE_3_bug_Au10": "Au 10 % по COMPOSITE_3 с ошибкой"}
    styles = {"Au10": ("#ff7f0e", "-"), "Ag10": ("#7f7f7f", "-"), "Cu10": ("#8c564b", "-"),
              "Au5+Ag5": ("#1f77b4", "--"), "Au5+Cu5": ("#d62728", "--"), "COMPOSITE_3_bug_Au10": ("k", ":")}
    for key, lab in labels.items():
        if ("jc", key) in extra:
            d = extra[("jc", key)]
            ax.plot(d[:, 0], d[:, 1], color=styles[key][0], ls=styles[key][1], lw=2, label=lab)
    ax.set_title("Два сорта наночастиц (Максвелл-Гарнетт, Джонсон-Кристи)", fontsize=9.5)
    ax.text(600, 0.62, "кривая с ошибкой уходит в P < 0\n(усиление вместо потерь), обрезано", fontsize=7.5, ha="center")
    ax.set_xlabel("длина волны в свободном пространстве, нм")
    ax.set_ylabel("поглощение P = 1 − R − T")
    ax.set_ylim(-0.02, 1.02)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8.5, loc="lower left")
    fig.suptitle("Решётка цилиндров композита Au 10 % (период 500 нм, D = 300 нм, h = 100 нм, слой 950 нм, "
                 "прослойка 325 нм, подложка n = 1,45), нормальное падение", fontsize=10.5)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    fig.savefig(L.RESULTS / "mixing_rules_on_grating.png", dpi=170)
    plt.close(fig)

    # схема ячейки: вид сбоку и вид сверху
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 4.2), gridspec_kw=dict(width_ratios=[1.6, 1]))
    z = 0
    a1.add_patch(Rectangle((-250, -1500), 500, 1500, color="#dfe7f0", ec="k", lw=0.8))
    a1.text(0, -750, "подложка n = 1,45\n(полубесконечная)", ha="center", va="center", fontsize=9)
    a1.add_patch(Rectangle((-250, 0), 500, 325, color="#eef3f8", ec="k", lw=0.8))
    a1.text(0, 162, "прослойка n = 1,45, 325 нм", ha="center", va="center", fontsize=9)
    a1.add_patch(Rectangle((-250, 325), 500, 950, color="#c9a36b", ec="k", lw=0.8))
    a1.text(0, 800, "композит Au 10 % в матрице n = 1,77, 950 нм", ha="center", va="center", fontsize=9)
    a1.add_patch(Rectangle((-150, 1275), 300, 100, color="#c9a36b", ec="k", lw=0.8))
    a1.text(0, 1325, "цилиндр D = 300", ha="center", va="center", fontsize=8.5)
    a1.text(-240, 1470, "воздух", ha="left", fontsize=9)
    a1.annotate("", xy=(0, 1420), xytext=(0, 1700), arrowprops=dict(arrowstyle="->", lw=1.5))
    a1.text(25, 1580, "плоская волна,\nнормальное падение", ha="left", va="center", fontsize=8.5)
    a1.annotate("", xy=(250, 1820), xytext=(-250, 1820), arrowprops=dict(arrowstyle="<->", lw=1))
    a1.text(0, 1850, "период 500 нм", ha="center", fontsize=8.5)
    a1.set_xlim(-300, 300); a1.set_ylim(-1500, 1950)
    a1.set_xlabel("x поперёк, нм"); a1.set_ylabel("z по нормали, нм")
    a1.set_title("Вид сбоку: одна ячейка решётки (сечение)")
    a1.set_aspect("auto")
    a2.add_patch(Rectangle((-250, -250), 500, 500, color="#f4f4f4", ec="k", lw=0.8))
    a2.add_patch(Ellipse((0, 0), 300, 300, color="#c9a36b", ec="k", lw=0.8))
    a2.text(0, 0, "композит\nD = 300 нм", ha="center", va="center", fontsize=9)
    a2.set_xlim(-300, 300); a2.set_ylim(-300, 300)
    a2.set_aspect("equal")
    a2.set_xlabel("x, нм"); a2.set_ylabel("y, нм")
    a2.set_title("Вид сверху: ячейка 500 × 500 нм")
    fig.suptitle("Метаповерхность из входного файла письма 19.08.2026: геометрия, которая решалась", fontsize=10.5)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(L.RESULTS / "mixing_rules_scheme.png", dpi=170)
    plt.close(fig)
    print("рисунки записаны")


if __name__ == "__main__":
    main()
