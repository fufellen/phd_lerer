from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import repeat_lerer_metal_strip as base


WIDTHS_NM = np.array([300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 1200.0, 1500.0, 2000.0])
LAMBDA_SWEEP_NM = np.arange(450.0, 801.0, 25.0)
WIDTH_PLOT_LAMBDAS_NM = [500.0, 600.0, 700.0]


def solve_all() -> tuple[dict[float, complex], dict[float, dict[float, complex]]]:
    planar = base.solve_planar_branch(LAMBDA_SWEEP_NM)
    by_width: dict[float, dict[float, complex]] = {}
    for width_nm in WIDTHS_NM:
        by_width[float(width_nm)] = base.solve_strip_branch(planar, width_nm / 1000.0)
    return planar, by_width


def write_csv(path: Path, planar: dict[float, complex], by_width: dict[float, dict[float, complex]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["lambda_nm", "width_nm", "neff_real", "neff_imag", "ratio_imag_real"])
        for lam, neff in planar.items():
            writer.writerow(
                [
                    f"{lam:.0f}",
                    "infinity",
                    f"{neff.real:.9g}",
                    f"{neff.imag:.9g}",
                    f"{base.attenuation_ratio(neff):.9g}",
                ]
            )
        for width_nm, branch in by_width.items():
            for lam, neff in branch.items():
                writer.writerow(
                    [
                        f"{lam:.0f}",
                        f"{width_nm:.0f}",
                        f"{neff.real:.9g}",
                        f"{neff.imag:.9g}",
                        f"{base.attenuation_ratio(neff):.9g}",
                    ]
                )


def plot_vs_lambda(path: Path, planar: dict[float, complex], by_width: dict[float, dict[float, complex]]) -> None:
    selected_widths = [400.0, 800.0, 1200.0, 2000.0]
    lambdas = np.array(list(planar.keys()))
    fig, (ax_loss, ax_neff) = plt.subplots(
        2, 1, figsize=(7.0, 6.5), sharex=True, gridspec_kw={"height_ratios": [1, 1.45]}
    )

    for width_nm in selected_widths:
        branch = by_width[width_nm]
        values = np.array([branch[float(lam)] for lam in lambdas])
        ax_loss.plot(lambdas, [base.attenuation_ratio(x) for x in values], marker="o", ms=3, label=f"W={width_nm:.0f} nm")
        ax_neff.plot(lambdas, values.real, marker="o", ms=3, label=f"W={width_nm:.0f} nm")

    inf_values = np.array([planar[float(lam)] for lam in lambdas])
    ax_loss.plot(lambdas, [base.attenuation_ratio(x) for x in inf_values], "k--", lw=1.5, label="W=infinity")
    ax_neff.plot(lambdas, inf_values.real, "k--", lw=1.5, label="W=infinity")

    ax_loss.set_ylabel("n''/n'")
    ax_loss.grid(True, color="0.85")
    ax_loss.legend(frameon=False, fontsize=8)
    ax_neff.set_ylabel("n'")
    ax_neff.set_xlabel("lambda, nm")
    ax_neff.grid(True, color="0.85")
    ax_neff.legend(frameon=False, fontsize=8)
    fig.suptitle("Analytic EDP-like Ag strip width sweep")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_vs_width(path: Path, planar: dict[float, complex], by_width: dict[float, dict[float, complex]]) -> None:
    widths = np.array(list(by_width.keys()))
    fig, (ax_neff, ax_loss) = plt.subplots(2, 1, figsize=(7.0, 6.0), sharex=True)

    for lam in WIDTH_PLOT_LAMBDAS_NM:
        values = np.array([by_width[float(width)][lam] for width in widths])
        ax_neff.plot(widths, values.real, marker="o", ms=3, label=f"lambda={lam:.0f} nm")
        ax_loss.plot(widths, [base.attenuation_ratio(x) for x in values], marker="o", ms=3, label=f"lambda={lam:.0f} nm")
        ax_neff.axhline(planar[lam].real, color="0.75", lw=0.8, ls="--")
        ax_loss.axhline(base.attenuation_ratio(planar[lam]), color="0.75", lw=0.8, ls="--")

    ax_neff.set_ylabel("n'")
    ax_neff.grid(True, color="0.85")
    ax_neff.legend(frameon=False, fontsize=8)
    ax_loss.set_ylabel("n''/n'")
    ax_loss.set_xlabel("W, nm")
    ax_loss.grid(True, color="0.85")
    ax_loss.legend(frameon=False, fontsize=8)
    fig.suptitle("Ag strip: dependence on finite width")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_summary(path: Path, planar: dict[float, complex], by_width: dict[float, dict[float, complex]]) -> None:
    lines = [
        "Analytic EDP-like width sweep for Ag strip on n=1.45 substrate",
        "Geometry: air | 20 nm Ag strip | substrate n=1.45.",
        "The W=infinity reference is the planar Ag film branch.",
        "",
        "lambda=500 nm:",
        "width_nm  n_real     n_imag     ratio",
    ]
    for width_nm in WIDTHS_NM:
        n = by_width[float(width_nm)][500.0]
        lines.append(f"{width_nm:8.0f}  {n.real:8.6f}  {n.imag:8.6f}  {base.attenuation_ratio(n):8.6f}")
    n_inf = planar[500.0]
    lines.append(f"{'inf':>8}  {n_inf.real:8.6f}  {n_inf.imag:8.6f}  {base.attenuation_ratio(n_inf):8.6f}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    out_dir = Path(__file__).resolve().parents[1] / "results"
    out_dir.mkdir(parents=True, exist_ok=True)
    planar, by_width = solve_all()
    write_csv(out_dir / "analytic_width_sweep_ag_strip.csv", planar, by_width)
    write_summary(out_dir / "analytic_width_sweep_ag_strip_summary.txt", planar, by_width)
    plot_vs_lambda(out_dir / "analytic_width_sweep_vs_lambda.png", planar, by_width)
    plot_vs_width(out_dir / "analytic_width_dependence_500_600_700nm.png", planar, by_width)


if __name__ == "__main__":
    main()
