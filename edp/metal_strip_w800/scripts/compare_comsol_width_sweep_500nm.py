from __future__ import annotations

import csv
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


LAMBDA_UM = 0.500
K0_UM = 2.0 * np.pi / LAMBDA_UM
WIDTHS_NM = [400, 800, 1200, 2000]


def parse_complex(text: str) -> complex:
    return complex(text.replace("i", "j"))


def load_analytic(path: Path) -> dict[int, complex]:
    out: dict[int, complex] = {}
    with path.open(encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if row["lambda_nm"] != "500" or row["width_nm"] == "infinity":
                continue
            width_nm = int(float(row["width_nm"]))
            out[width_nm] = complex(float(row["neff_real"]), float(row["neff_imag"]))
    return out


def load_comsol_modes(path: Path) -> list[complex]:
    modes: list[complex] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("%") or not re.match(r"^[0-9.+-]", line):
            continue
        cols = next(csv.reader([line]))
        try:
            modes.append(parse_complex(cols[0]))
        except ValueError:
            continue
    return modes


def attenuation_ratio(neff: complex) -> float:
    return abs(neff.imag / neff.real)


def power_length_um(neff: complex) -> float:
    return 1.0 / (2.0 * K0_UM * abs(neff.imag))


def select_mode(modes: list[complex], target: complex) -> tuple[complex, int]:
    candidates = [m for m in modes if m.real > 1.45 and abs(m.imag) > 1e-3]
    pool = candidates if candidates else modes
    selected = min(pool, key=lambda m: abs(m - target))
    return selected, len(candidates)


def main() -> None:
    project_dir = Path(__file__).resolve().parents[1]
    results_dir = project_dir / "results"
    comsol_dir = project_dir / "comsol_exports"
    results_dir.mkdir(parents=True, exist_ok=True)
    analytic = load_analytic(results_dir / "analytic_width_sweep_ag_strip.csv")
    rows = []
    for width_nm in WIDTHS_NM:
        target = analytic[width_nm]
        modes = load_comsol_modes(comsol_dir / f"comsol_ag_strip_width_sweep_500nm_w{width_nm:04d}_modes.csv")
        selected, candidate_count = select_mode(modes, target)
        rows.append(
            {
                "width_nm": width_nm,
                "analytic_real": target.real,
                "analytic_imag": target.imag,
                "analytic_ratio": attenuation_ratio(target),
                "analytic_L_power_um": power_length_um(target),
                "comsol_real": selected.real,
                "comsol_imag": selected.imag,
                "comsol_ratio": attenuation_ratio(selected),
                "comsol_L_power_um": power_length_um(selected),
                "delta_real": selected.real - target.real,
                "delta_ratio": attenuation_ratio(selected) - attenuation_ratio(target),
                "mode_count": len(modes),
                "candidate_count": candidate_count,
            }
        )

    fields = [
        "width_nm",
        "analytic_real",
        "analytic_imag",
        "analytic_ratio",
        "analytic_L_power_um",
        "comsol_real",
        "comsol_imag",
        "comsol_ratio",
        "comsol_L_power_um",
        "delta_real",
        "delta_ratio",
        "mode_count",
        "candidate_count",
    ]
    with (results_dir / "comsol_width_sweep_500nm_selected_modes.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: (f"{v:.12g}" if isinstance(v, float) else v) for k, v in row.items()})

    lines = [
        "COMSOL width sweep selected modes at lambda=500 nm",
        "Selection rule: Re(neff)>1.45 and |Im(neff)|>1e-3, then closest complex neff to analytic EDP-like target.",
        "",
        "width_nm  analytic_neff              comsol_neff                analytic_ratio  comsol_ratio  analytic_Lum  comsol_Lum  delta_Re_percent  delta_ratio_percent",
    ]
    for row in rows:
        delta_re_percent = row["delta_real"] / row["analytic_real"] * 100.0
        delta_ratio_percent = row["delta_ratio"] / row["analytic_ratio"] * 100.0
        lines.append(
            f"{row['width_nm']:8d}  "
            f"{row['analytic_real']:.9f}+{row['analytic_imag']:.9f}i  "
            f"{row['comsol_real']:.9f}+{row['comsol_imag']:.9f}i  "
            f"{row['analytic_ratio']:.9f}  {row['comsol_ratio']:.9f}  "
            f"{row['analytic_L_power_um']:.4f}  {row['comsol_L_power_um']:.4f}  "
            f"{delta_re_percent:16.3f}  {delta_ratio_percent:19.3f}"
        )
    (results_dir / "comsol_width_sweep_500nm_selected_modes_summary.txt").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )

    widths = np.array([row["width_nm"] for row in rows], dtype=float)
    fig, axes = plt.subplots(3, 1, figsize=(7.0, 8.0), sharex=True)
    axes[0].plot(widths, [row["analytic_real"] for row in rows], "o-", label="EDP-like analytic")
    axes[0].plot(widths, [row["comsol_real"] for row in rows], "s-", label="COMSOL selected mode")
    axes[0].set_ylabel("n'")
    axes[0].grid(True, color="0.85")
    axes[0].legend(frameon=False, fontsize=8)

    axes[1].plot(widths, [row["analytic_ratio"] for row in rows], "o-", label="EDP-like analytic")
    axes[1].plot(widths, [row["comsol_ratio"] for row in rows], "s-", label="COMSOL selected mode")
    axes[1].set_ylabel("n''/n'")
    axes[1].grid(True, color="0.85")
    axes[1].legend(frameon=False, fontsize=8)

    axes[2].plot(widths, [row["analytic_L_power_um"] for row in rows], "o-", label="EDP-like analytic")
    axes[2].plot(widths, [row["comsol_L_power_um"] for row in rows], "s-", label="COMSOL selected mode")
    axes[2].set_ylabel("L_power, um")
    axes[2].set_xlabel("W, nm")
    axes[2].grid(True, color="0.85")
    axes[2].legend(frameon=False, fontsize=8)

    fig.suptitle("Ag strip width sweep at lambda=500 nm")
    fig.tight_layout()
    fig.savefig(results_dir / "comsol_vs_analytic_width_sweep_500nm.png", dpi=180)
    plt.close(fig)

    print("\n".join(lines))


if __name__ == "__main__":
    main()
