"""Print the numeric facts quoted in the manuscript text.

Reads the same CSVs as make_figures.py, so text numbers always match the
plotted curves. Usage: python article_os/scripts/metrics_for_text.py
"""

import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from make_figures import read_csv, xy


def peak(lams, vals, lo=None, hi=None):
    pts = [
        (l, v)
        for l, v in zip(lams, vals)
        if (lo is None or l >= lo) and (hi is None or l <= hi)
    ]
    l, v = max(pts, key=lambda p: p[1])
    return l, v


def uv_mean(lams, vals, cut=400.0):
    pts = [v for l, v in zip(lams, vals) if l <= cut]
    return sum(pts) / len(pts)


def half_max_band(lams, vals):
    vmax = max(vals)
    half = vmax / 2.0
    inside = [l for l, v in zip(lams, vals) if v >= half]
    return min(inside), max(inside), max(inside) - min(inside)


def frohlich(lams, re_eps, eps_m):
    """Wavelength where Re eps_p + 2*eps_m crosses zero (linear interp)."""
    f = [r + 2 * eps_m for r in re_eps]
    for i in range(1, len(f)):
        if f[i - 1] * f[i] < 0:
            t = f[i - 1] / (f[i - 1] - f[i])
            return lams[i - 1] + t * (lams[i] - lams[i - 1])
    return None


d1 = read_csv("grating_2d_rcwa/results/reproduce_fig2_au_vs_cu.csv")
d2 = read_csv("grating_2d_rcwa/results/reproduce_fig2_au_vs_cu_2d.csv")
d5 = read_csv("grating_2d_rcwa/results/reproduce_fig5_hfn_vs_tin.csv")
dr = read_csv("vibr_t/results/rerun_fig2_lerer_constants.csv")

print("== 1D strips, J&C constants (fig. 4 solid) ==")
for mat, lc, pc in (("Au", "lambda_nm", "P_Au"), ("Cu", "lambda_nm_Cu", "P_Cu")):
    lam, p = xy(d1, lc, pc)
    pk_vis = peak(lam, p, lo=450, hi=700)
    pk_all = peak(lam, p)
    print(
        f"{mat}: visible peak {pk_vis[0]:.0f} nm P={pk_vis[1]:.3f}; "
        f"global peak {pk_all[0]:.0f} nm P={pk_all[1]:.3f}; "
        f"UV mean(<=400) {uv_mean(lam, p):.3f}"
    )

print("== 2D cylinders, J&C constants (fig. 4 dashed) ==")
for mat, lc, pc in (("Au", "lambda_nm", "P_Au"), ("Cu", "lambda_nm_Cu", "P_Cu")):
    lam, p = xy(d2, lc, pc)
    pk_vis = peak(lam, p, lo=450, hi=700)
    pk_all = peak(lam, p)
    print(
        f"{mat}: visible peak {pk_vis[0]:.0f} nm P={pk_vis[1]:.3f}; "
        f"global peak {pk_all[0]:.0f} nm P={pk_all[1]:.3f}; "
        f"UV mean(<=400, from {min(lam):.0f}) {uv_mean(lam, p):.3f}"
    )

print("== HfN vs TiN, 1D (fig. 5) ==")
for mat, lc, pc in (("HfN", "lambda_nm_HfN", "P_HfN"), ("TiN", "lambda_nm_TiN", "P_TiN")):
    lam, p = xy(d5, lc, pc)
    pk = peak(lam, p)
    lo, hi, width = half_max_band(lam, p)
    print(
        f"{mat}: peak {pk[0]:.0f} nm P={pk[1]:.3f}; "
        f"band P>=half-max: {lo:.0f}-{hi:.0f} nm, width {width:.0f} nm"
    )

print("== two constants bases, 1D (fig. 7) ==")
for col in ("P_Au_lerer", "P_Cu_lerer", "P_Au_jc", "P_Cu_jc"):
    lam, p = xy(dr, "lambda_nm", col)
    pk_vis = peak(lam, p, lo=450, hi=700)
    print(
        f"{col}: UV mean(<=400) {uv_mean(lam, p):.3f}; "
        f"visible peak {pk_vis[0]:.0f} nm P={pk_vis[1]:.3f}"
    )

print("== Frohlich crossings (host n=1.77) ==")
dc = read_csv("vibr_t/results/compare_constants_au_cu.csv")
eps_m = 1.77 * 1.77
for col, name in (
    ("Au_re_lerer", "Au Lerer"),
    ("Au_re_jc", "Au J&C"),
    ("Cu_re_lerer", "Cu Lerer (extrapolated <517)"),
    ("Cu_re_jc", "Cu J&C"),
):
    lam, re = xy(dc, "lambda_nm", col)
    lf = frohlich(lam, re, eps_m)
    print(f"{name}: Frohlich at {lf:.1f} nm" if lf else f"{name}: no crossing")
