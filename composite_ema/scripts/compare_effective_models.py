"""Compare effective-medium variants for Au/Cu nanoparticle composites.

This script is a diagnostic companion to ``composite_ema.py``.  It keeps the
same Lerer sign convention as the original COMPOSIT code:

    eps = eps' - i eps'',  Im(eps) <= 0 for passive media.

Models compared:
- Maxwell-Garnett / Clausius-Mossotti with bulk particle permittivity.
- Bruggeman / CPA for the same host and inclusion.
- Maxwell-Garnett through MLWA-corrected dipole polarizability for finite
  spherical particles with selected radii.

The MLWA expression is evaluated in the common optics convention
``exp(-i omega t)`` and converted back by complex conjugation.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

import composite_ema as ce


RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"


@dataclass(frozen=True)
class ModelResult:
    material: str
    model: str
    lambda_nm: float
    eps_eff: complex


def bruggeman_from_eps(eps_h: complex, eps_i: complex, f: float) -> complex:
    """Symmetric Bruggeman EMA for two spherical components.

    Equation:

        f (eps_i - eps)/(eps_i + 2 eps)
      + (1-f) (eps_h - eps)/(eps_h + 2 eps) = 0.

    The physical branch is chosen as the root closest to the host permittivity.
    For small f this is the root that continuously tends to eps_h.
    """
    b = f * (2.0 * eps_i - eps_h) + (1.0 - f) * (2.0 * eps_h - eps_i)
    roots = np.roots([2.0, -b, -eps_i * eps_h])
    return min(roots, key=lambda z: abs(z - eps_h))


def mlwa_mg_from_eps(
    eps_h_lerer: complex,
    eps_p_lerer: complex,
    f: float,
    lam_nm: float,
    radius_nm: float,
) -> complex:
    """MG through MLWA-corrected polarizability of a finite sphere.

    The current material tables use the Lerer/SVF convention
    ``eps = eps' - i eps''``.  The MLWA formula below is written for
    ``exp(-i omega t)`` where passive loss has positive Im(eps), so the
    calculation is done on conjugated permittivities and then conjugated back.
    """
    eps_h = eps_h_lerer.conjugate()
    eps_p = eps_p_lerer.conjugate()

    radius_m = radius_nm * 1e-9
    lam_m = lam_nm * 1e-9
    n_h = complex(np.sqrt(eps_h)).real
    k_h = 2.0 * math.pi * n_h / lam_m

    beta0 = radius_m**3 * (eps_p - eps_h) / (eps_p + 2.0 * eps_h)
    beta_mlwa = beta0 / (
        1.0
        - (k_h**2 / radius_m) * beta0
        - 1j * (2.0 / 3.0) * (k_h**3) * beta0
    )
    eta = f * beta_mlwa / radius_m**3
    eps_eff = eps_h * (1.0 + 2.0 * eta) / (1.0 - eta)
    return eps_eff.conjugate()


def material_eps(material: str, lam_nm: float) -> complex:
    if material == "Au":
        return ce.eps_au(lam_nm)
    if material == "Cu":
        return ce.eps_cu(lam_nm)
    raise ValueError(material)


def compute_rows() -> list[ModelResult]:
    n_host = 1.77
    eps_h = complex(n_host**2, 0.0)
    f = 0.10
    lam_grid = np.arange(300.0, 900.0 + 1e-9, 2.0)
    radii_nm = (10.0, 30.0, 50.0)

    rows: list[ModelResult] = []
    for material in ("Au", "Cu"):
        for lam_nm in lam_grid:
            eps_p = material_eps(material, float(lam_nm))
            rows.append(
                ModelResult(material, "MG bulk / COMPOSIT", float(lam_nm), ce.composite_from_eps(eps_h, eps_p, f))
            )
            rows.append(
                ModelResult(material, "Bruggeman / CPA", float(lam_nm), bruggeman_from_eps(eps_h, eps_p, f))
            )
            for radius_nm in radii_nm:
                rows.append(
                    ModelResult(
                        material,
                        f"MG + MLWA sphere R={radius_nm:.0f} nm",
                        float(lam_nm),
                        mlwa_mg_from_eps(eps_h, eps_p, f, float(lam_nm), radius_nm),
                    )
                )
    return rows


def summarize(rows: list[ModelResult]) -> list[dict[str, object]]:
    summary: list[dict[str, object]] = []
    keys = sorted({(row.material, row.model) for row in rows})
    for material, model in keys:
        subset = [row for row in rows if row.material == material and row.model == model]
        peak = max(subset, key=lambda row: abs(row.eps_eff.imag))
        sign_violations = sum(1 for row in subset if row.eps_eff.imag > 1e-9)
        summary.append(
            {
                "material": material,
                "model": model,
                "peak_lambda_nm": peak.lambda_nm,
                "peak_Re_eps_eff": peak.eps_eff.real,
                "peak_Im_eps_eff": peak.eps_eff.imag,
                "max_abs_Im_eps_eff": abs(peak.eps_eff.imag),
                "positive_Im_points": sign_violations,
            }
        )
    return summary


def write_outputs(rows: list[ModelResult], summary: list[dict[str, object]]) -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    csv_path = RESULTS_DIR / "effective_model_comparison_AuCu.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["material", "model", "lambda_nm", "Re_eps_eff", "Im_eps_eff", "loss_value_minus_Im"],
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "material": row.material,
                    "model": row.model,
                    "lambda_nm": f"{row.lambda_nm:.1f}",
                    "Re_eps_eff": f"{row.eps_eff.real:.8g}",
                    "Im_eps_eff": f"{row.eps_eff.imag:.8g}",
                    "loss_value_minus_Im": f"{-row.eps_eff.imag:.8g}",
                }
            )

    summary_path = RESULTS_DIR / "effective_model_comparison_AuCu_summary.csv"
    with summary_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(summary[0].keys()))
        writer.writeheader()
        for row in summary:
            writer.writerow(row)

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
    colors = {
        "MG bulk / COMPOSIT": "#1f77b4",
        "Bruggeman / CPA": "#ff7f0e",
        "MG + MLWA sphere R=10 nm": "#2ca02c",
        "MG + MLWA sphere R=30 nm": "#9467bd",
        "MG + MLWA sphere R=50 nm": "#d62728",
    }
    for ax, material in zip(axes, ("Au", "Cu")):
        for model in colors:
            subset = [row for row in rows if row.material == material and row.model == model]
            ax.plot(
                [row.lambda_nm for row in subset],
                [-row.eps_eff.imag for row in subset],
                label=model,
                color=colors[model],
                linewidth=1.2,
            )
        ax.set_ylabel(f"{material}: -Im eps_eff")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=7, ncol=2)
    axes[-1].set_xlabel("wavelength, nm")
    fig.suptitle("Effective-medium variants, n_host=1.77, f=10%, Lerer sign convention", fontsize=11)
    fig.tight_layout()
    png_path = RESULTS_DIR / "effective_model_comparison_AuCu.png"
    fig.savefig(png_path, dpi=160)

    print(f"[write] {csv_path}")
    print(f"[write] {summary_path}")
    print(f"[write] {png_path}")


def main() -> None:
    rows = compute_rows()
    summary = summarize(rows)
    write_outputs(rows, summary)

    print("[summary]")
    for row in summary:
        print(
            f"{row['material']:>2} | {row['model']:<28} | "
            f"peak {row['peak_lambda_nm']:>5.1f} nm | "
            f"Im={row['peak_Im_eps_eff']:.4g} | "
            f"positive-Im points={row['positive_Im_points']}"
        )


if __name__ == "__main__":
    main()
