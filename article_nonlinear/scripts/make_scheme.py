"""
Схема исследуемой метаповерхности для нелинейной статьи
(«Интенсивностно-зависимое поглощение метаповерхностей с ПНЧ»).

Геометрия та же, что в линейной работе (article_os/scripts/make_figures.py,
fig1_scheme), но схема несёт нелинейный смысл: падающая волна помечена
интенсивностью I, композитный слой - интенсивностно-зависимой
eps_eff(lambda, I), а врезка под стопкой показывает золотые наночастицы
с eps_p(I) = eps_L + Delta eps_NL.

Подписи на английском - как в фигурах линейной статьи.

Запуск (из корня репозитория):
    python article_nonlinear/scripts/make_scheme.py
Выход: article_nonlinear/figures/fig_scheme.{png,eps}
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle

BLACK = "black"
OUT = Path(__file__).resolve().parents[1] / "figures"


def panel_letter(ax, letter: str, x: float = 0.0, y: float = 1.0) -> None:
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=12,
            fontstyle="italic", va="top")


def main() -> None:
    fig, (ax_top, ax_side) = plt.subplots(
        1, 2, figsize=(7.6, 3.6), layout="constrained",
        gridspec_kw=dict(width_ratios=[1.0, 1.55]),
    )

    # ---- (a) вид сверху: квадратная решётка цилиндров -------------------
    ax_top.set_aspect("equal")
    for i in range(3):
        for j in range(3):
            ax_top.add_patch(Circle((i, j), 0.4, fill=False, color=BLACK, lw=1.1))
    # период d
    ax_top.annotate("", xy=(2.0, -0.70), xytext=(1.0, -0.70),
                    arrowprops=dict(arrowstyle="<->", color=BLACK, lw=0.9))
    ax_top.text(1.5, -0.98, "d", ha="center", fontstyle="italic", fontsize=10)
    # диаметр D (внутри верхнего левого цилиндра)
    ax_top.annotate("", xy=(0.4, 2.0), xytext=(-0.4, 2.0),
                    arrowprops=dict(arrowstyle="<->", color=BLACK, lw=0.9))
    ax_top.text(0.0, 2.14, "D", ha="center", fontstyle="italic", fontsize=10)
    ax_top.set_xlim(-0.85, 2.75)
    ax_top.set_ylim(-1.25, 2.55)
    ax_top.axis("off")
    panel_letter(ax_top, "a")

    # ---- (b) поперечное сечение стопки + врезка снизу --------------------
    ax_side.set_aspect("equal")
    x0, w = 0.0, 3.0
    ax_side.add_patch(Rectangle((x0, 0.0), w, 0.55, facecolor="0.85", edgecolor=BLACK, lw=0.9))
    ax_side.add_patch(Rectangle((x0, 0.55), w, 0.35, facecolor="0.62", edgecolor=BLACK, lw=0.9))
    for k in range(4):
        ax_side.add_patch(Rectangle((x0 + 0.12 + 0.78 * k, 0.90), 0.55, 0.35,
                                    facecolor="0.62", edgecolor=BLACK, lw=0.9))

    # падающая волна с интенсивностью I (в промежутке между столбиками)
    ax_side.annotate("", xy=(x0 + 1.56, 1.38), xytext=(x0 + 1.56, 2.02),
                     arrowprops=dict(arrowstyle="->", color=BLACK, lw=1.2))
    ax_side.text(x0 + 1.70, 2.06, "incident wave, intensity $I$", fontsize=8.5, va="center")

    # подписи справа с выносками
    lx = x0 + w + 0.30
    for y_src, y_txt, text in [
        (1.45, 1.72, "air"),
        (1.07, 1.20, "pillars (composite)"),
        (0.72, 0.68, r"composite, $\varepsilon_{\rm eff}(\lambda, I)$"),
        (0.27, 0.16, "substrate"),
    ]:
        ax_side.annotate(text, xy=(x0 + w, y_src), xytext=(lx, y_txt),
                         va="center", fontsize=8.5,
                         arrowprops=dict(arrowstyle="-", color=BLACK, lw=0.7))

    # врезка снизу: золотые наночастицы в матрице.
    # Выноска начинается ВНУТРИ композитного слоя (маркер) и рисуется поверх
    # подложки (zorder), иначе читается, будто увеличена подложка.
    cx, cy, r = x0 + 1.5, -0.95, 0.40
    y_src_inset = 0.72  # середина композитного слоя (0.55..0.90)
    ax_side.plot([cx, cx], [y_src_inset, cy + r], color="0.35", lw=0.8, ls=":", zorder=5)
    ax_side.plot([cx], [y_src_inset], marker="o", ms=3.0, color="0.15", zorder=6)
    ax_side.add_patch(Circle((cx, cy), r, facecolor="0.62", edgecolor=BLACK, lw=0.9))
    for dx, dy in [(-0.17, 0.11), (0.09, 0.19), (0.17, -0.09), (-0.09, -0.17), (0.0, 0.0)]:
        ax_side.add_patch(Circle((cx + dx, cy + dy), 0.052, facecolor=BLACK, edgecolor="none"))
    ax_side.text(cx + r + 0.18, cy + 0.10, r"Au NPs, $C=10\%$", fontsize=8.5, va="center")
    ax_side.text(cx + r + 0.18, cy - 0.18,
                 r"$\varepsilon_p(I)=\varepsilon_L+\Delta\varepsilon_{\rm NL}$",
                 fontsize=8.5, va="center")

    ax_side.set_xlim(-0.25, 5.85)
    ax_side.set_ylim(-1.65, 2.35)
    ax_side.axis("off")
    panel_letter(ax_side, "b")

    OUT.mkdir(parents=True, exist_ok=True)
    for ext, kw in (("png", dict(dpi=600)), ("eps", {})):
        p = OUT / f"fig_scheme.{ext}"
        fig.savefig(p, **kw)
        print(f"записан: {p}")


if __name__ == "__main__":
    main()
