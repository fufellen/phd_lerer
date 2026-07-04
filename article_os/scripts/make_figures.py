"""Generate journal-ready figures for the "Optika i spektroskopiya" manuscript.

Rules implemented here (journals.ioffe.ru author rules):
- axis labels, units and all in-figure text in English only;
- curves are numbered (1, 2, ...) and decoded in the caption, not in a legend;
- no grid inside plots; grayscale-friendly line styles (solid/dashed/dotted);
- each figure saved as a separate file, PNG (600 dpi) + EPS (vector).

Inputs: CSV results already present in this repository (composite_ema,
grating_2d_rcwa, vibr_t). The script recomputes nothing.

Usage:
    python article_os/scripts/make_figures.py
Outputs land in article_os/figures/.
"""

import csv
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
OUT = os.path.join(REPO, "article_os", "figures")

BLACK = "black"
GRAY = "0.45"
STYLES = ["-", "--", ":", "-.", (0, (5, 1, 1, 1, 1, 1))]

plt.rcParams.update(
    {
        "font.size": 11,
        "axes.linewidth": 0.8,
        "lines.linewidth": 1.4,
        "savefig.bbox": "tight",
    }
)


def read_csv(rel_path):
    """Read a CSV into {column: [floats / strings / None for empty cells]}.

    Columns may be ragged (e.g. the Au and Cu wavelength grids in the 2D
    run start at different wavelengths), so values are converted one by one.
    """
    path = os.path.join(REPO, rel_path)
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    cols = {}
    for key in rows[0]:
        vals = []
        for r in rows:
            v = r[key]
            if v is None or v.strip() == "":
                vals.append(None)
                continue
            try:
                vals.append(float(v))
            except ValueError:
                vals.append(v)
        cols[key] = vals
    return cols


def xy(d, xcol, ycol):
    """Paired columns with empty cells dropped."""
    pairs = [
        (x, y) for x, y in zip(d[xcol], d[ycol]) if x is not None and y is not None
    ]
    return [p[0] for p in pairs], [p[1] for p in pairs]


def save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    for ext in ("png", "eps"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"), dpi=600)
    plt.close(fig)
    print(f"saved {name}.png / {name}.eps")


def label_curve(ax, x, y, text, dx=0, dy=0):
    """Put a curve number near the point (x, y)."""
    ax.annotate(
        text,
        xy=(x, y),
        xytext=(x + dx, y + dy),
        fontsize=11,
        fontstyle="italic",
    )


def panel_letter(ax, letter):
    ax.text(
        0.02,
        0.97,
        letter,
        transform=ax.transAxes,
        fontsize=12,
        fontstyle="italic",
        va="top",
    )


# ---------------------------------------------------------------- fig. 1
def fig1_scheme():
    fig, (ax_top, ax_side) = plt.subplots(1, 2, figsize=(7.2, 3.2), layout="constrained")

    # (a) top view: square lattice of cylinders, period 500 nm, D = 400 nm
    ax_top.set_aspect("equal")
    for i in range(3):
        for j in range(3):
            ax_top.add_patch(
                Circle((i * 1.0, j * 1.0), 0.4, fill=False, color=BLACK, lw=1.2)
            )
    # period marker
    ax_top.annotate(
        "",
        xy=(2.0, -0.62),
        xytext=(1.0, -0.62),
        arrowprops=dict(arrowstyle="<->", color=BLACK, lw=0.9),
    )
    ax_top.text(1.5, -0.85, "d", ha="center", fontstyle="italic")
    ax_top.annotate(
        "",
        xy=(0.4, 2.0),
        xytext=(-0.4, 2.0),
        arrowprops=dict(arrowstyle="<->", color=BLACK, lw=0.9),
    )
    ax_top.text(0.0, 2.18, "D", ha="center", fontstyle="italic")
    ax_top.set_xlim(-0.8, 2.8)
    ax_top.set_ylim(-1.05, 2.6)
    ax_top.axis("off")
    panel_letter(ax_top, "a")

    # (b) cross-section of the stack
    ax_side.set_aspect("equal")
    x0, w = 0.0, 3.0
    # substrate
    ax_side.add_patch(
        Rectangle((x0, 0.0), w, 0.55, fill=True, facecolor="0.85", edgecolor=BLACK)
    )
    # continuous composite layer
    ax_side.add_patch(
        Rectangle((x0, 0.55), w, 0.35, fill=True, facecolor="0.65", edgecolor=BLACK)
    )
    # grating pillars (same composite)
    for k in range(4):
        ax_side.add_patch(
            Rectangle(
                (x0 + 0.12 + 0.78 * k, 0.90),
                0.55,
                0.35,
                fill=True,
                facecolor="0.65",
                edgecolor=BLACK,
            )
        )
    labels = [
        (0.27, 0.05, "substrate"),
        (0.72, 0.45, "composite layer"),
        (1.07, 0.95, "composite pillars"),
        (1.60, 1.60, "air"),
    ]
    for y_src, y_txt, text in labels:
        ax_side.annotate(
            text,
            xy=(x0 + w, y_src),
            xytext=(x0 + w + 0.55, y_txt),
            va="center",
            fontsize=10,
            arrowprops=dict(arrowstyle="-", color=BLACK, lw=0.7),
        )
    ax_side.annotate(
        "",
        xy=(x0 + 1.5, 1.40),
        xytext=(x0 + 1.5, 1.95),
        arrowprops=dict(arrowstyle="->", color=BLACK, lw=1.1),
    )
    ax_side.text(x0 + 1.7, 2.05, "incident wave", fontsize=10)
    ax_side.set_xlim(-0.2, 6.4)
    ax_side.set_ylim(-0.15, 2.3)
    ax_side.axis("off")
    panel_letter(ax_side, "b")

    save(fig, "fig1")


# ---------------------------------------------------------------- fig. 2
def fig2_eps_eff():
    d = read_csv("composite_ema/results/composite_eps_Au_vs_Cu.csv")
    lam = d["lambda_nm"]
    fig, (ax_re, ax_im) = plt.subplots(1, 2, figsize=(7.2, 3.1), layout="constrained")

    ax_re.plot(lam, d["Re_eps_eff_Au_C10"], linestyle=STYLES[0], color=BLACK)
    ax_re.plot(lam, d["Re_eps_eff_Cu_C10"], linestyle=STYLES[1], color=BLACK)
    ax_re.set_xlabel("$\\lambda$, nm")
    ax_re.set_ylabel("Re $\\varepsilon_{eff}$")
    ax_re.set_xlim(250, 1000)
    label_curve(ax_re, 560, 7.7, "1")
    label_curve(ax_re, 660, 4.6, "2")
    panel_letter(ax_re, "a")

    ax_im.plot(lam, d["Im_eps_eff_Au_C10"], linestyle=STYLES[0], color=BLACK)
    ax_im.plot(lam, d["Im_eps_eff_Cu_C10"], linestyle=STYLES[1], color=BLACK)
    ax_im.set_xlabel("$\\lambda$, nm")
    ax_im.set_ylabel("Im $\\varepsilon_{eff}$")
    ax_im.set_xlim(250, 1000)
    label_curve(ax_im, 545, -5.7, "1")
    label_curve(ax_im, 640, -2.8, "2")
    panel_letter(ax_im, "b")

    save(fig, "fig2")


# ---------------------------------------------------------------- fig. 3
def fig3_ema_models():
    d = read_csv("composite_ema/results/effective_model_comparison_AuCu.csv")
    models = [
        ("MG bulk", lambda m: m.startswith("MG bulk")),
        ("Bruggeman", lambda m: m.startswith("Bruggeman")),
        ("MLWA R=10", lambda m: "MLWA" in m and "R=10" in m),
        ("MLWA R=30", lambda m: "MLWA" in m and "R=30" in m),
        ("MLWA R=50", lambda m: "MLWA" in m and "R=50" in m),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.1), layout="constrained")
    for ax, mat, letter in ((axes[0], "Au", "a"), (axes[1], "Cu", "b")):
        for idx, (_, match) in enumerate(models):
            lam = [
                d["lambda_nm"][i]
                for i in range(len(d["material"]))
                if d["material"][i] == mat and match(d["model"][i])
            ]
            im = [
                d["Im_eps_eff"][i]
                for i in range(len(d["material"]))
                if d["material"][i] == mat and match(d["model"][i])
            ]
            if not lam:
                raise ValueError(f"no rows matched for {mat} model index {idx}")
            color = BLACK if idx < 2 else GRAY
            ax.plot(lam, im, linestyle=STYLES[idx % len(STYLES)], color=color)
        ax.set_xlabel("$\\lambda$, nm")
        ax.set_ylabel("Im $\\varepsilon_{eff}$")
        panel_letter(ax, letter)
    # curve numbers, placed by peak positions from the summary CSV
    label_curve(axes[0], 540, -6.3, "1")
    label_curve(axes[0], 735, -2.9, "2")
    label_curve(axes[0], 600, -6.9, "3")
    label_curve(axes[0], 665, -6.3, "4")
    label_curve(axes[0], 815, -3.6, "5")
    label_curve(axes[1], 565, -3.5, "1")
    label_curve(axes[1], 730, -2.7, "2")
    label_curve(axes[1], 615, -4.0, "3")
    label_curve(axes[1], 665, -5.7, "4")
    label_curve(axes[1], 810, -3.4, "5")
    save(fig, "fig3")


# ---------------------------------------------------------------- fig. 4
def fig4_p_1d_vs_2d():
    d1 = read_csv("grating_2d_rcwa/results/reproduce_fig2_au_vs_cu.csv")
    d2 = read_csv("grating_2d_rcwa/results/reproduce_fig2_au_vs_cu_2d.csv")
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.1), layout="constrained")

    ax = axes[0]
    ax.plot(*xy(d1, "lambda_nm", "P_Au"), linestyle=STYLES[0], color=BLACK)
    ax.plot(*xy(d2, "lambda_nm", "P_Au"), linestyle=STYLES[1], color=BLACK)
    ax.set_xlabel("$\\lambda$, nm")
    ax.set_ylabel("P")
    ax.set_ylim(0, 1)
    label_curve(ax, 430, 0.66, "1")
    label_curve(ax, 560, 0.95, "2")
    panel_letter(ax, "a")

    ax = axes[1]
    ax.plot(*xy(d1, "lambda_nm_Cu", "P_Cu"), linestyle=STYLES[0], color=BLACK)
    ax.plot(*xy(d2, "lambda_nm_Cu", "P_Cu"), linestyle=STYLES[1], color=BLACK)
    ax.set_xlabel("$\\lambda$, nm")
    ax.set_ylabel("P")
    ax.set_ylim(0, 1)
    label_curve(ax, 380, 0.86, "1")
    label_curve(ax, 470, 0.70, "2")
    panel_letter(ax, "b")

    save(fig, "fig4")


# ---------------------------------------------------------------- fig. 5
def fig5_hfn_tin():
    d = read_csv("grating_2d_rcwa/results/reproduce_fig5_hfn_vs_tin.csv")
    fig, ax = plt.subplots(figsize=(4.2, 3.2))
    ax.plot(d["lambda_nm_HfN"], d["P_HfN"], linestyle=STYLES[0], color=BLACK)
    ax.plot(d["lambda_nm_TiN"], d["P_TiN"], linestyle=STYLES[1], color=BLACK)
    ax.set_xlabel("$\\lambda$, nm")
    ax.set_ylabel("P")
    ax.set_ylim(0, 1)
    label_curve(ax, 610, 0.82, "1")
    label_curve(ax, 800, 0.60, "2")
    save(fig, "fig5")


# ---------------------------------------------------------------- fig. 6
def fig6_constants():
    d = read_csv("vibr_t/results/compare_constants_au_cu.csv")
    lam = d["lambda_nm"]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.1), layout="constrained")

    ax = axes[0]
    ax.plot(lam, d["Au_re_lerer"], linestyle=STYLES[0], color=BLACK)
    ax.plot(lam, d["Au_re_jc"], linestyle=STYLES[1], color=BLACK)
    ax.set_xlabel("$\\lambda$, nm")
    ax.set_ylabel("Re $\\varepsilon$")
    ax.set_xlim(200, 900)
    ax.set_ylim(-35, 15)
    label_curve(ax, 340, 1.5, "1")
    label_curve(ax, 640, -13.5, "2")
    panel_letter(ax, "a")

    ax = axes[1]
    # shade the Cu extrapolation region (no tabulated data below 517 nm)
    ax.axvspan(200, 517, facecolor="0.88", edgecolor="none")
    ax.plot(lam, d["Cu_re_lerer"], linestyle=STYLES[0], color=BLACK)
    ax.plot(lam, d["Cu_re_jc"], linestyle=STYLES[1], color=BLACK)
    ax.set_xlabel("$\\lambda$, nm")
    ax.set_ylabel("Re $\\varepsilon$")
    ax.set_xlim(200, 900)
    ax.set_ylim(-35, 15)
    label_curve(ax, 330, 6.5, "1")
    label_curve(ax, 350, -4.8, "2")
    panel_letter(ax, "b")

    save(fig, "fig6")


# ---------------------------------------------------------------- fig. 7
def fig7_p_two_bases():
    d = read_csv("vibr_t/results/rerun_fig2_lerer_constants.csv")
    lam = d["lambda_nm"]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.1), layout="constrained")

    ax = axes[0]
    ax.plot(lam, d["P_Au_lerer"], linestyle=STYLES[0], color=BLACK)
    ax.plot(lam, d["P_Cu_lerer"], linestyle=STYLES[1], color=BLACK)
    ax.set_xlabel("$\\lambda$, nm")
    ax.set_ylabel("P")
    ax.set_ylim(0, 1)
    label_curve(ax, 250, 0.88, "1")
    label_curve(ax, 330, 0.55, "2")
    panel_letter(ax, "a")

    ax = axes[1]
    ax.plot(lam, d["P_Au_jc"], linestyle=STYLES[0], color=BLACK)
    ax.plot(lam, d["P_Cu_jc"], linestyle=STYLES[1], color=BLACK)
    ax.set_xlabel("$\\lambda$, nm")
    ax.set_ylabel("P")
    ax.set_ylim(0, 1)
    label_curve(ax, 520, 0.92, "1")
    label_curve(ax, 400, 0.84, "2")
    panel_letter(ax, "b")

    save(fig, "fig7")


if __name__ == "__main__":
    fig1_scheme()
    fig2_eps_eff()
    fig3_ema_models()
    fig4_p_1d_vs_2d()
    fig5_hfn_tin()
    fig6_constants()
    fig7_p_two_bases()
    print("done")
