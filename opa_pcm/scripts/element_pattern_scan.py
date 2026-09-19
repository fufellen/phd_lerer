"""Pattern of a real emitter in the 1D OPA: scan loss and sidelobes at N = 32.

The array model of simulate_pcm_opa.py is an array factor of isotropic elements. Here the power
pattern is multiplied by the pattern of a real emitter - a silicon wire 450 x 220 nm ending in air,
pitch 0.775 um = lambda/2 at 1.55 um - computed by the full-wave COMSOL model of the emitter array
(RunArrayF0, 19.09.2026; plane-wave spectrum of the field 0.3 um in front of the facets, both
transverse components):
  * isolated  - the same wire with no neighbours;
  * embedded  - the mean of the two middle channels of a four-emitter array, each solved with its
                own port excited and the neighbours passive (channel 2 is the mirror image of
                channel 1, so the mean is symmetric).
The multiplication is the usual large-array approximation: every element is taken to radiate with
the same pattern. Ideal shifters only: the element pattern multiplies the power pattern whatever the
shifter does, so its effect on the metrics of the other scenarios is the same factor.

Printed and written to results/element_pattern_scan_n32.csv for target angles -60...60 deg:
pointing error, highest sidelobe (outside the main lobe as in simulate_pcm_opa.pattern_metrics),
and the scan loss - the power in the target direction relative to the broadside beam of the same
array. Figure: results/element_pattern_scan_n32.png.

Run: python opa_pcm/scripts/element_pattern_scan.py   (NumPy and Matplotlib)
"""

from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import simulate_pcm_opa as S  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "element_patterns_z300.csv"
OUT = ROOT / "results"
N = 32
TARGETS = np.arange(-60.0, 60.0001, 5.0)


def load_patterns() -> dict[str, np.ndarray]:
    raw = np.loadtxt(DATA, delimiter=",", skiprows=1)
    theta = raw[:, 0]

    def on_grid(col: np.ndarray) -> np.ndarray:
        # outside the computed range (|theta| > 89 deg) the emitter is taken not to radiate
        return np.interp(S.THETA_DEG, theta, col, left=0.0, right=0.0)

    iso = on_grid(raw[:, 1])
    emb = on_grid(0.5 * (raw[:, 3] + raw[:, 4]))
    i0 = int(np.argmin(np.abs(S.THETA_DEG)))
    return {"isotropic": np.ones_like(S.THETA_DEG),
            "isolated": iso / iso[i0],
            "embedded": emb / emb[i0]}


def sweep(ep: np.ndarray) -> list[dict]:
    rows = []
    broadside = None
    for t0 in TARGETS:
        q = S.target_phases(N, float(t0))
        p = S.array_power(q, np.ones(N)) * ep
        m = S.pattern_metrics(p, np.ones(N), float(t0), N)
        if broadside is None:
            q0 = S.target_phases(N, 0.0)
            broadside = float(np.max(S.array_power(q0, np.ones(N)) * ep))
        target = float(np.interp(t0, S.THETA_DEG, p))
        rows.append({"target_deg": float(t0), "peak_deg": m.peak_angle_deg,
                     "pointing_error_deg": m.pointing_error_deg, "sll_db": m.sll_db,
                     "scan_loss_db": 10.0 * math.log10(max(target, 1e-300) / broadside)})
    return rows


def main() -> None:
    eps = load_patterns()
    for name in ("isolated", "embedded"):
        e = eps[name]
        above = S.THETA_DEG[e >= 0.5]
        print("%-9s element: -3 dB width %.1f deg, level at 30 / 45 / 60 deg: %s dB" % (
            name, above.max() - above.min(),
            " / ".join("%.1f" % (10 * math.log10(max(float(np.interp(a, S.THETA_DEG, e)), 1e-12)))
                       for a in (30.0, 45.0, 60.0))))
    results = {name: sweep(ep) for name, ep in eps.items()}
    OUT.mkdir(exist_ok=True)
    with open(OUT / "element_pattern_scan_n32.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["element", "target_deg", "peak_deg", "pointing_error_deg", "sll_db", "scan_loss_db"])
        for name, rows in results.items():
            for r in rows:
                w.writerow([name, r["target_deg"], "%.3f" % r["peak_deg"], "%.3f" % r["pointing_error_deg"],
                            "%.2f" % r["sll_db"], "%.2f" % r["scan_loss_db"]])
    print("\n%7s | %-26s | %-26s | %-26s" % ("target", "isotropic: err, SLL, loss", "isolated",
                                             "embedded"))
    for k, t0 in enumerate(TARGETS):
        if t0 < 0:
            continue
        print("%6.0f° | %s" % (t0, " | ".join(
            "%+6.2f° %6.1f %6.2f dB" % (results[n][k]["pointing_error_deg"], results[n][k]["sll_db"],
                                         results[n][k]["scan_loss_db"])
            for n in ("isotropic", "isolated", "embedded"))))

    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
    styles = {"isotropic": ("k-", "1 — изотропный элемент"),
              "isolated": ("b--", "2 — одиночная жила"),
              "embedded": ("r-.", "3 — жила в решётке")}
    for name, ep in eps.items():
        st, lab = styles[name]
        ax[0].plot(S.THETA_DEG, 10 * np.log10(np.maximum(ep, 1e-6)), st, lw=1.3, label=lab)
        rows = results[name]
        ax[1].plot(TARGETS, [r["scan_loss_db"] for r in rows], st, lw=1.3, label=lab)
        ax[2].plot(TARGETS, [r["sll_db"] for r in rows], st, lw=1.3, label=lab)
    ax[0].set(xlim=(-90, 90), ylim=(-20, 1), xlabel="угол, град", ylabel="диаграмма элемента, дБ")
    ax[1].set(xlabel="заданный угол, град", ylabel="потери сканирования, дБ")
    ax[2].set(xlabel="заданный угол, град", ylabel="наибольший боковой лепесток, дБ")
    for a in ax:
        a.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "element_pattern_scan_n32.png", dpi=150)
    print("\nwritten:", OUT / "element_pattern_scan_n32.csv", OUT / "element_pattern_scan_n32.png")


if __name__ == "__main__":
    main()
