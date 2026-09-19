"""Journal figures for the nonlinear metasurface manuscript («Оптика и спектроскопия»).

Same rules as article_os/scripts/make_figures.py (journals.ioffe.ru author rules):
- axis labels, units and all in-figure text in English only;
- curves are numbered (1, 2, ...) and decoded in the caption, not in a legend;
- no grid inside plots; grayscale-friendly line styles;
- each figure saved as a separate file, PNG (600 dpi) + EPS (vector).

The script recomputes nothing except the Frohlich wavelength marked in fig2 (the same helper
as nonlinear_au_sweep.py). Inputs are CSV results already in this repository:
  fig2     eps_eff(lambda) of the composite at four intensities
           composite_ema/results/nonlinear_au_sweep_spectra.csv
  fig3     absorption of the grating at several intensities, 1D approximation, decoupled
           grating_2d_rcwa/results/rcwa_nonlinear_au_absorption.csv
  fig3_sc  2D grating of cylinders: decoupled vs self-consistent absorption shift and the
           mean field in the composite (variant of fig3 for the self-consistent version)
           metasurface_grcwa/results/nonlinear_selfconsistent.csv
  fig4     two-temperature dynamics for a 2 ps pulse
           composite_ema/results/ttm_hot_electron_dynamics.csv
  fig5     peak electron and lattice temperatures vs pulse duration
           composite_ema/results/ttm_hot_electron_regimes.csv

Usage (from the repository root):
    python article_nonlinear/scripts/make_figures.py
Outputs land in article_nonlinear/figures/. fig1 (the scheme) is made by make_scheme.py.
"""

import csv
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
OUT = os.path.join(REPO, "article_nonlinear", "figures")

BLACK = "black"
GRAY = "0.45"
INTENSITIES = (0.0, 1e6, 5e6, 2e7)  # W/cm^2, the series of the manuscript

plt.rcParams.update(
    {
        "font.size": 11,
        "axes.linewidth": 0.8,
        "lines.linewidth": 1.4,
        "savefig.bbox": "tight",
    }
)


def rows(rel_path):
    with open(os.path.join(REPO, rel_path), newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    for ext in ("png", "eps"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"), dpi=600)
    plt.close(fig)
    print(f"saved {name}.png / {name}.eps")


def label_curve(ax, x, y, text, xy=None, ha="left"):
    """Curve number at (x, y); with xy, a thin leader line to the curve point xy."""
    if xy is None:
        ax.annotate(text, xy=(x, y), fontsize=11, fontstyle="italic", ha=ha)
    else:
        ax.annotate(
            text,
            xy=xy,
            xytext=(x, y),
            fontsize=11,
            fontstyle="italic",
            ha=ha,
            arrowprops=dict(arrowstyle="-", lw=0.6, color=BLACK, shrinkA=1, shrinkB=1),
        )


def panel_letter(ax, letter, x=0.02, y=0.97):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=12, fontstyle="italic", va="top")


def frohlich_nm():
    """Frohlich condition Re eps_Au = -2 eps_m for the host of the manuscript, as in fig. 2 of the
    working note (nonlinear_au_sweep.py)."""
    sys.path.insert(0, os.path.join(REPO, "composite_ema", "scripts"))
    import nonlinear_au_sweep as nl

    return nl.find_frohlich_crossing(nl.eps_au, nl.N_HOST ** 2, 400.0, 900.0).lam_nm


# ---------------------------------------------------------------- fig. 2
def fig2_eps_eff():
    # I = 1e6 W/cm^2 is left out: on this scale it coincides with the linear curve
    data = rows("composite_ema/results/nonlinear_au_sweep_spectra.csv")
    lam = [float(r["lambda_nm"]) for r in data]
    keys = ["0e+00", "5e+06", "2e+07"]
    styles = [dict(color=BLACK, ls="-"), dict(color=BLACK, ls="--"),
              dict(color=BLACK, ls=":", lw=1.8)]
    lam_f = frohlich_nm()
    fig, (ax_re, ax_im) = plt.subplots(2, 1, figsize=(4.4, 5.4), sharex=True, layout="constrained")
    for k, st in zip(keys, styles):
        ax_re.plot(lam, [float(r["Re_eps_eff_I" + k]) for r in data], **st)
        ax_im.plot(lam, [float(r["Im_eps_eff_I" + k]) for r in data], **st)
    for ax in (ax_re, ax_im):
        ax.axvline(lam_f, color=GRAY, ls=(0, (1, 3)), lw=0.8)
        ax.set_xlim(400, 900)
    ax_re.set_ylim(0.9, 8.2)
    ax_re.set_ylabel("Re $\\varepsilon_{eff}$")
    ax_im.set_ylabel("Im $\\varepsilon_{eff}$")
    ax_im.set_xlabel("$\\lambda$, nm")
    # Re: the curves part at the dip near 546 nm (1.74 / 1.95 / 2.31); 1 from below,
    # 3 from above, 2 from above across the dotted curve 3
    label_curve(ax_re, 575, 1.05, "1", xy=(549, 1.77))
    label_curve(ax_re, 492, 3.35, "2", xy=(541, 1.93))
    label_curve(ax_re, 532, 3.55, "3", xy=(547, 2.33))
    # Im: at the minimum near 574 nm (-6.04 / -5.52 / -4.65)
    label_curve(ax_im, 640, -6.25, "1", xy=(580, -5.6))
    label_curve(ax_im, 640, -4.95, "2", xy=(583, -5.05))
    label_curve(ax_im, 640, -3.85, "3", xy=(586, -4.35))
    panel_letter(ax_re, "a", x=0.93)
    panel_letter(ax_im, "b", x=0.93, y=0.30)
    save(fig, "fig2")
    return lam_f


# ---------------------------------------------------------------- fig. 3
def fig3_p_1d():
    data = rows("grating_2d_rcwa/results/rcwa_nonlinear_au_absorption.csv")
    tab = {}
    for r in data:
        tab.setdefault(float(r["I_wcm2"]), []).append((float(r["lambda_nm"]), float(r["P"])))
    for v in tab.values():
        v.sort()
    lam = [p[0] for p in tab[0.0]]
    p0 = [p[1] for p in tab[0.0]]
    fig, (ax_p, ax_d) = plt.subplots(2, 1, figsize=(4.4, 5.4), sharex=True, layout="constrained")
    ax_p.plot(lam, p0, color=BLACK, ls="-")
    ax_p.plot(lam, [p[1] for p in tab[2e7]], color=BLACK, ls="--")
    ax_p.set_ylabel("$P$")
    ax_p.set_ylim(0, 1)
    # the two spectra part only at 545-600 nm: 1 below, 2 above
    label_curve(ax_p, 552, 0.64, "1", xy=(576, 0.766))
    label_curve(ax_p, 590, 0.93, "2", xy=(572, 0.806))
    ax_d.axhline(0, color=GRAY, lw=0.6)
    styles = {1e6: dict(color=BLACK, ls=":", lw=1.8), 5e6: dict(color=BLACK, ls="--"),
              2e7: dict(color=BLACK, ls="-")}
    for i_w, st in styles.items():
        ax_d.plot(lam, [p[1] - q for p, q in zip(tab[i_w], p0)], **st)
    ax_d.set_ylabel("$\\Delta P$")
    ax_d.set_xlabel("$\\lambda$, nm")
    ax_d.set_xlim(400, 900)
    ax_d.set_ylim(-0.013, 0.042)
    # numbers stacked over the maxima at 565 nm (0.0039 / 0.0153 / 0.0353)
    label_curve(ax_d, 566, 0.0050, "1", ha="center")
    label_curve(ax_d, 566, 0.0167, "2", ha="center")
    label_curve(ax_d, 566, 0.0368, "3", ha="center")
    panel_letter(ax_p, "a", x=0.93)
    panel_letter(ax_d, "b", x=0.93)
    save(fig, "fig3")


# ---------------------------------------------------------- fig. 3, self-consistent variant
def fig3_selfconsistent():
    data = rows("metasurface_grcwa/results/nonlinear_selfconsistent.csv")
    top = [r for r in data if float(r["I_wcm2"]) == 2e7]
    lam = [float(r["lambda_nm"]) for r in top]
    col = lambda key: [float(r[key]) for r in top]
    fig, (ax_p, ax_d, ax_e) = plt.subplots(3, 1, figsize=(4.4, 7.2), sharex=True, layout="constrained")
    ax_p.plot(lam, col("P_lin"), color=BLACK)
    ax_p.set_ylabel("$P$")
    ax_p.set_ylim(0, 1)
    ax_d.axhline(0, color=GRAY, lw=0.6)
    ax_d.plot(lam, col("dP_dec"), color=BLACK, ls="--")
    ax_d.plot(lam, col("dP_sc"), color=BLACK, ls="-")
    ax_d.set_ylabel("$\\Delta P$")
    # 1 left of the rise of the dashed curve to its maximum at 610 nm, 2 over the maximum of the
    # solid one at 560 nm
    label_curve(ax_d, 592, 0.0235, "1", ha="center")
    label_curve(ax_d, 560, 0.0135, "2", ha="center")
    ax_e.axhline(1, color=GRAY, ls=":", lw=1.0)
    ax_e.plot(lam, col("ratio_cyl_lin"), color=BLACK, ls="-")
    ax_e.plot(lam, col("ratio_film_lin"), color=BLACK, ls="--")
    ax_e.set_ylabel("$n_h\\langle|E|^2\\rangle/|E_0|^2$")
    ax_e.set_xlabel("$\\lambda$, nm")
    ax_e.set_xlim(400, 700)
    ax_e.set_ylim(0, 1.3)
    label_curve(ax_e, 520, 0.50, "1", xy=(545, 0.64))
    label_curve(ax_e, 470, 0.22, "2", xy=(500, 0.42))
    panel_letter(ax_p, "a")
    panel_letter(ax_d, "b")
    panel_letter(ax_e, "c")
    save(fig, "fig3_sc")


# ---------------------------------------------------------------- fig. 4
def fig4_ttm_dynamics():
    data = rows("composite_ema/results/ttm_hot_electron_dynamics.csv")
    t = [float(r["t_ps"]) for r in data]
    fig, ax = plt.subplots(figsize=(4.4, 3.2), layout="constrained")
    ax.plot(t, [float(r["T_e_K"]) for r in data], color=BLACK, ls="-")
    ax.plot(t, [float(r["T_l_K"]) for r in data], color=BLACK, ls="--")
    ax.set_xlim(0, 20)
    ax.set_xlabel("$t$, ps")
    ax.set_ylabel("$T$, K")
    label_curve(ax, 7.2, 560, "1", xy=(5.9, 515))
    # 2 just over the lattice curve where the electron curve is still far above it
    label_curve(ax, 6.3, 318, "2", ha="center")
    save(fig, "fig4")


# ---------------------------------------------------------------- fig. 5
def fig5_ttm_regimes():
    data = rows("composite_ema/results/ttm_hot_electron_regimes.csv")
    tau = [float(r["tau_fwhm_ps"]) for r in data]
    fig, ax = plt.subplots(figsize=(4.4, 3.2), layout="constrained")
    ax.semilogx(tau, [float(r["Te_peak_K"]) for r in data], color=BLACK, ls="-", marker="o", ms=4)
    ax.semilogx(tau, [float(r["Tl_peak_K"]) for r in data], color=BLACK, ls="--", marker="s", ms=4,
                mfc="white")
    ax.axvline(1.0, color=GRAY, ls=":", lw=1.0)
    ax.set_xlabel("$\\tau_p$, ps")
    ax.set_ylabel("$T_{max}$, K")
    label_curve(ax, 12, 520, "1", xy=(6.5, 455))
    label_curve(ax, 0.3, 330, "2", ha="center")
    save(fig, "fig5")


if __name__ == "__main__":
    lam_f = fig2_eps_eff()
    print(f"Frohlich wavelength marked in fig2: {lam_f:.1f} nm")
    fig3_p_1d()
    fig3_selfconsistent()
    fig4_ttm_dynamics()
    fig5_ttm_regimes()
