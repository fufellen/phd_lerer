"""Рисунки к pillar_shapes.py: схема форм столбиков, спектры и средние по полосе.
Выход: results/pillar_shapes_scheme.png, results/pillar_shapes.png."""
from __future__ import annotations

import csv
import io
import math
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Rectangle, Ellipse, Polygon  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import lerer_grcwa as L  # noqa: E402
import pillar_shapes as PS  # noqa: E402

COL = {"square": "#1f77b4", "circle_area": "#2ca02c", "circle_side": "#98df8a", "ellipse_x": "#9467bd",
       "ellipse_y": "#c5b0d5", "frustum": "#ff7f0e", "cone": "#d62728", "cone_circle": "#e377c2",
       "cone_h150": "#8c564b", "frustum_h150": "#c49c94", "two_step": "#17becf"}
SHORT = {"square": "квадрат", "circle_area": "круг той же площади", "circle_side": "круг D = сторона",
         "ellipse_x": "эллипс, E ∥ большой оси", "ellipse_y": "эллипс, E ∥ малой оси", "frustum": "усечённый конус",
         "cone": "конус", "cone_circle": "круглый конус", "cone_h150": "конус, h = 150", "frustum_h150": "усеч. конус, h = 150",
         "two_step": "две ступени"}


def read(path, col):
    rows = list(csv.DictReader(io.open(path, encoding="utf-8")))
    out = {}
    for r in rows:
        out.setdefault((r["design"], r["variant"]), []).append((float(r["lambda_nm"]), float(r[col])))
    return {k: np.array(v) for k, v in out.items()}


def scheme():
    fig, axs = plt.subplots(2, 4, figsize=(12, 5.6))
    a = 125.0
    tops = [("квадрат 125 × 125", lambda ax: ax.add_patch(Rectangle((-62.5, -62.5), 125, 125, color="#c9a36b", ec="k"))),
            ("круг той же площади, D = 141", lambda ax: ax.add_patch(Ellipse((0, 0), 141, 141, color="#c9a36b", ec="k"))),
            ("эллипс той же площади 173 × 115", lambda ax: ax.add_patch(Ellipse((0, 0), 172.7, 115.1, color="#c9a36b", ec="k"))),
            ("круг D = 125 (меньше площадь)", lambda ax: ax.add_patch(Ellipse((0, 0), 125, 125, color="#c9a36b", ec="k")))]
    for ax, (t, draw) in zip(axs[0], tops):
        ax.add_patch(Rectangle((-100, -100), 200, 200, color="#f4f4f4", ec="k", lw=0.8))
        draw(ax)
        ax.set_xlim(-110, 110); ax.set_ylim(-110, 110); ax.set_aspect("equal")
        ax.set_title("вид сверху: " + t, fontsize=9)
        ax.set_xlabel("x, нм"); ax.set_ylabel("y, нм")
    sides = [("цилиндр (столбик) 125 × 100", [(-62.5, 0), (62.5, 0), (62.5, 100), (-62.5, 100)]),
             ("усечённый конус 125 → 62", [(-62.5, 0), (62.5, 0), (31, 100), (-31, 100)]),
             ("конус 125 → 12", [(-62.5, 0), (62.5, 0), (6, 100), (-6, 100)]),
             ("две ступени 125 и 75, по 100 нм", [(-62.5, 0), (62.5, 0), (62.5, 100), (37.5, 100), (37.5, 200), (-37.5, 200), (-37.5, 100), (-62.5, 100)])]
    for ax, (t, pts) in zip(axs[1], sides):
        ax.add_patch(Rectangle((-100, -60), 200, 60, color="#c9a36b", ec="k", lw=0.8))
        ax.text(0, -30, "композит", ha="center", va="center", fontsize=8)
        ax.add_patch(Polygon(pts, closed=True, color="#c9a36b", ec="k", lw=0.8))
        if "конус" in t:
            n = 8
            for i in range(1, n):
                zc = 100 * i / n
                ax.axhline(zc, xmin=0.3, xmax=0.7, color="k", lw=0.4, ls=":")
        ax.set_xlim(-110, 110); ax.set_ylim(-60, 220); ax.set_aspect("equal")
        ax.set_title("вид сбоку: " + t, fontsize=9)
        ax.set_xlabel("x поперёк, нм"); ax.set_ylabel("z по нормали, нм")
    fig.suptitle("Формы столбиков для предложенной конструкции (период 200 нм); для исходной все размеры вдвое больше. "
                 "Конус решается лестницей из 8 слоёв (пунктир)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(L.RESULTS / "pillar_shapes_scheme.png", dpi=170)
    plt.close(fig)


def main():
    scheme()
    fine = read(L.RESULTS / "pillar_shapes_fine.csv", "A_layers")
    series = read(L.RESULTS / "pillar_shapes_series.csv", "A_layers")
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.8), gridspec_kw=dict(width_ratios=[1.2, 1.2, 1.4]))
    ax = axs[0]
    for v in ("square", "circle_area", "frustum", "cone", "two_step", "ellipse_x", "ellipse_y"):
        if ("improved", v) in fine:
            d = fine[("improved", v)]
            ax.plot(d[:, 0], d[:, 1], color=COL[v], lw=1.8, ls="--" if v == "ellipse_y" else "-", label=SHORT[v])
    ax.set_title("Предложенная конструкция: период 200 нм,\nдва слоя композита, золотое зеркало", fontsize=10)
    ax.set_ylim(0.95, 1.001)
    ax = axs[1]
    for v in ("square", "cone", "two_step"):
        if ("original", v) in fine:
            d = fine[("original", v)]
            ax.plot(d[:, 0], d[:, 1], color=COL[v], lw=1.8, label=SHORT[v])
    ax.set_title("Исходная конструкция: период 400 нм,\nкомпозит 2000 нм, подложка n = 1,45", fontsize=10)
    ax.set_ylim(0.5, 1.001)
    for ax in axs[:2]:
        ax.set_xlabel("длина волны в свободном пространстве, нм")
        ax.set_ylabel("поглощение в слоях, 1 − R − T")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8.5, loc="lower left")
    ax = axs[2]
    variants = PS.VARIANTS
    x = np.arange(len(variants))
    for k, (design, off, col) in enumerate((("improved", -0.2, "#1f77b4"), ("original", 0.2, "#ff7f0e"))):
        means = [np.mean(series[(design, v)][:, 1]) if (design, v) in series else np.nan for v in variants]
        ax.bar(x + off, means, width=0.4, color=col, label="предложенная (период 200)" if design == "improved" else "исходная (период 400)")
        for xi, m in zip(x + off, means):
            if not np.isnan(m):
                ax.text(xi, m + 0.002, "%.3f" % m, ha="center", fontsize=6.5, rotation=90)
    ax.set_xticks(x)
    ax.set_xticklabels([SHORT[v] for v in variants], rotation=60, ha="right", fontsize=8)
    ax.set_ylim(0.85, 1.03)
    ax.set_ylabel("среднее поглощение по 8 длинам волн 400…750 нм")
    ax.set_title("Средние по полосе для всех форм", fontsize=10)
    ax.grid(axis="y", alpha=0.3)
    ax.legend(fontsize=8.5, loc="upper right")
    fig.tight_layout()
    fig.savefig(L.RESULTS / "pillar_shapes.png", dpi=170)
    plt.close(fig)
    print("рисунки записаны")


if __name__ == "__main__":
    main()
