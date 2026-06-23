from __future__ import annotations

import csv
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


LAMBDAS_NM = [450, 500, 550, 600, 650, 700, 750, 800]


def parse_complex(text: str) -> complex:
    return complex(text.replace("i", "j"))


def attenuation_ratio(neff: complex) -> float:
    return abs(neff.imag / neff.real)


def power_length_um(lambda_nm: int, neff: complex) -> float:
    lambda_um = lambda_nm / 1000.0
    return lambda_um / (4.0 * np.pi * abs(neff.imag))


def load_analytic(path: Path) -> dict[int, complex]:
    out: dict[int, complex] = {}
    with path.open(encoding="utf-8") as f:
        for row in csv.DictReader(f):
            lambda_nm = int(float(row["lambda_nm"]))
            if lambda_nm in LAMBDAS_NM:
                out[lambda_nm] = complex(float(row["neff_w800_real"]), float(row["neff_w800_imag"]))
    return out


def load_comsol_modes(path: Path) -> list[dict[str, float | complex | int]]:
    modes: list[dict[str, float | complex | int]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("%") or not re.match(r"^[0-9.+-]", line):
            continue
        cols = next(csv.reader([line]))
        if len(cols) < 13:
            continue
        try:
            modes.append(
                {
                    "mode_row": len(modes) + 1,
                    "neff": parse_complex(cols[0]),
                    "neff_real": float(cols[2]),
                    "neff_imag": float(cols[3]),
                    "ratio": abs(float(cols[4])),
                    "beta_rad_um": float(cols[5]),
                    "damp_db_um_raw": float(cols[6]),
                    "field_ag": float(cols[7]),
                    "field_sub": float(cols[8]),
                    "field_air": float(cols[9]),
                    "ag_frac": float(cols[10]),
                    "sub_frac": float(cols[11]),
                    "air_frac": float(cols[12]),
                }
            )
        except ValueError:
            continue
    return modes


def choose_mode(
    lambda_nm: int,
    modes: list[dict[str, float | complex | int]],
    target: complex,
    previous: complex | None,
) -> tuple[dict[str, float | complex | int], str, int, int]:
    physical = [
        m
        for m in modes
        if float(m["neff_real"]) > 1.45 and abs(float(m["neff_imag"])) > 1e-3
    ]
    field_pool = [
        m
        for m in physical
        if np.isfinite(float(m["ag_frac"]))
        and float(m["ag_frac"]) >= 1e-2
        and 0.05 <= float(m["sub_frac"]) <= 0.95
    ]
    pool = field_pool or physical or modes
    rule = "field_pool_closest_to_edp"
    if not field_pool and physical:
        rule = "no_field_pool__physical_closest_to_edp"
    elif not physical:
        rule = "fallback_all_modes_closest_to_edp"

    # The field pool rejects light-line dielectric modes. The analytic EDP target then
    # picks the branch; the continuity term is deliberately weak and only breaks ties.
    lambda_scale = max(0.02, abs(target))

    def score(mode: dict[str, float | complex | int]) -> float:
        neff = mode["neff"]
        assert isinstance(neff, complex)
        edp_score = abs(neff - target) / lambda_scale
        continuity = 0.0
        if previous is not None:
            continuity = abs(neff - previous) / max(0.02, abs(previous))
        return edp_score + 0.15 * continuity

    selected = min(pool, key=score)
    return selected, rule, len(physical), len(field_pool)


def main() -> None:
    project_dir = Path(__file__).resolve().parents[1]
    results_dir = project_dir / "results"
    comsol_dir = project_dir / "comsol_exports"
    results_dir.mkdir(parents=True, exist_ok=True)
    analytic = load_analytic(results_dir / "analytic_metal_strip_w800.csv")
    selected_rows: list[dict[str, float | int | str]] = []
    candidate_rows: list[dict[str, float | int | str]] = []
    previous: complex | None = None

    for lambda_nm in LAMBDAS_NM:
        csv_path = comsol_dir / f"comsol_ag_strip_w800_lambda_sweep_lam{lambda_nm:04d}nm_modes.csv"
        if not csv_path.exists():
            continue
        target = analytic[lambda_nm]
        modes = load_comsol_modes(csv_path)
        if not modes:
            continue
        selected, rule, physical_count, field_pool_count = choose_mode(lambda_nm, modes, target, previous)
        neff = selected["neff"]
        assert isinstance(neff, complex)
        previous = neff

        for mode in modes:
            m_neff = mode["neff"]
            assert isinstance(m_neff, complex)
            if float(mode["neff_real"]) > 1.45 and abs(float(mode["neff_imag"])) > 1e-3:
                candidate_rows.append(
                    {
                        "lambda_nm": lambda_nm,
                        "mode_row": int(mode["mode_row"]),
                        "neff_real": m_neff.real,
                        "neff_imag": m_neff.imag,
                        "ratio": attenuation_ratio(m_neff),
                        "ag_frac": float(mode["ag_frac"]),
                        "sub_frac": float(mode["sub_frac"]),
                        "air_frac": float(mode["air_frac"]),
                        "distance_to_edp": abs(m_neff - target),
                    }
                )

        selected_rows.append(
            {
                "lambda_nm": lambda_nm,
                "analytic_real": target.real,
                "analytic_imag": target.imag,
                "analytic_ratio": attenuation_ratio(target),
                "analytic_L_power_um": power_length_um(lambda_nm, target),
                "comsol_real": neff.real,
                "comsol_imag": neff.imag,
                "comsol_ratio": attenuation_ratio(neff),
                "comsol_L_power_um": power_length_um(lambda_nm, neff),
                "delta_real": neff.real - target.real,
                "delta_imag": neff.imag - target.imag,
                "delta_ratio": attenuation_ratio(neff) - attenuation_ratio(target),
                "delta_real_percent": (neff.real - target.real) / target.real * 100.0,
                "delta_ratio_percent": (
                    (attenuation_ratio(neff) - attenuation_ratio(target))
                    / attenuation_ratio(target)
                    * 100.0
                ),
                "mode_row": int(selected["mode_row"]),
                "mode_count": len(modes),
                "physical_candidate_count": physical_count,
                "field_pool_count": field_pool_count,
                "ag_frac": float(selected["ag_frac"]),
                "sub_frac": float(selected["sub_frac"]),
                "air_frac": float(selected["air_frac"]),
                "selection_rule": rule,
            }
        )

    if not selected_rows:
        raise SystemExit("No COMSOL lambda-sweep CSV files were found.")

    selected_fields = [
        "lambda_nm",
        "analytic_real",
        "analytic_imag",
        "analytic_ratio",
        "analytic_L_power_um",
        "comsol_real",
        "comsol_imag",
        "comsol_ratio",
        "comsol_L_power_um",
        "delta_real",
        "delta_imag",
        "delta_ratio",
        "delta_real_percent",
        "delta_ratio_percent",
        "mode_row",
        "mode_count",
        "physical_candidate_count",
        "field_pool_count",
        "ag_frac",
        "sub_frac",
        "air_frac",
        "selection_rule",
    ]
    with (results_dir / "comsol_lambda_sweep_w800_selected_modes.csv").open(
        "w", newline="", encoding="utf-8"
    ) as f:
        writer = csv.DictWriter(f, fieldnames=selected_fields)
        writer.writeheader()
        for row in selected_rows:
            writer.writerow({k: (f"{v:.12g}" if isinstance(v, float) else v) for k, v in row.items()})

    candidate_fields = [
        "lambda_nm",
        "mode_row",
        "neff_real",
        "neff_imag",
        "ratio",
        "ag_frac",
        "sub_frac",
        "air_frac",
        "distance_to_edp",
    ]
    with (results_dir / "comsol_lambda_sweep_w800_candidate_modes.csv").open(
        "w", newline="", encoding="utf-8"
    ) as f:
        writer = csv.DictWriter(f, fieldnames=candidate_fields)
        writer.writeheader()
        for row in candidate_rows:
            writer.writerow({k: (f"{v:.12g}" if isinstance(v, float) else v) for k, v in row.items()})

    lines = [
        "COMSOL wavelength sweep selected modes for Ag strip W=800 nm",
        (
            "Selection rule: reject light-line modes by field proxy "
            "(Re(neff)>1.45, |Im(neff)|>1e-3, Ag fraction >= 1e-2, 0.05<=substrate fraction<=0.95), "
            "then select the closest complex neff to the analytic EDP-like target with a weak continuity tie-breaker."
        ),
        "Field fractions are diagnostic |E|^2 integrals, not normalized modal power.",
        "",
        "lambda_nm  analytic_neff              comsol_neff                ratio_edp  ratio_comsol  Lp_edp_um  Lp_comsol_um  delta_Re_%  delta_ratio_%  mode_row  ag/sub/air fractions",
    ]
    for row in selected_rows:
        lines.append(
            f"{int(row['lambda_nm']):9d}  "
            f"{float(row['analytic_real']):.9f}+{float(row['analytic_imag']):.9f}i  "
            f"{float(row['comsol_real']):.9f}+{float(row['comsol_imag']):.9f}i  "
            f"{float(row['analytic_ratio']):.9f}  {float(row['comsol_ratio']):.9f}  "
            f"{float(row['analytic_L_power_um']):.4f}  {float(row['comsol_L_power_um']):.4f}  "
            f"{float(row['delta_real_percent']):10.3f}  {float(row['delta_ratio_percent']):13.3f}  "
            f"{int(row['mode_row']):8d}  "
            f"{float(row['ag_frac']):.3f}/{float(row['sub_frac']):.3f}/{float(row['air_frac']):.3f}"
        )
    (results_dir / "comsol_lambda_sweep_w800_selected_modes_summary.txt").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )

    lambdas = np.array([float(row["lambda_nm"]) for row in selected_rows])
    fig, axes = plt.subplots(4, 1, figsize=(7.2, 9.2), sharex=True)
    axes[0].plot(lambdas, [float(row["analytic_real"]) for row in selected_rows], "o-", label="EDP-like")
    axes[0].plot(lambdas, [float(row["comsol_real"]) for row in selected_rows], "s-", label="COMSOL")
    axes[0].set_ylabel("n'")
    axes[0].grid(True, color="0.85")
    axes[0].legend(frameon=False, fontsize=8)

    axes[1].plot(lambdas, [float(row["analytic_ratio"]) for row in selected_rows], "o-", label="EDP-like")
    axes[1].plot(lambdas, [float(row["comsol_ratio"]) for row in selected_rows], "s-", label="COMSOL")
    axes[1].set_ylabel("n''/n'")
    axes[1].grid(True, color="0.85")
    axes[1].legend(frameon=False, fontsize=8)

    axes[2].plot(lambdas, [float(row["analytic_L_power_um"]) for row in selected_rows], "o-", label="EDP-like")
    axes[2].plot(lambdas, [float(row["comsol_L_power_um"]) for row in selected_rows], "s-", label="COMSOL")
    axes[2].set_ylabel("L_power, um")
    axes[2].grid(True, color="0.85")
    axes[2].legend(frameon=False, fontsize=8)

    axes[3].plot(lambdas, [float(row["ag_frac"]) for row in selected_rows], "o-", label="Ag")
    axes[3].plot(lambdas, [float(row["sub_frac"]) for row in selected_rows], "s-", label="substrate")
    axes[3].plot(lambdas, [float(row["air_frac"]) for row in selected_rows], "^-", label="air")
    axes[3].set_xlabel("lambda, nm")
    axes[3].set_ylabel("|E|^2 fraction proxy")
    axes[3].grid(True, color="0.85")
    axes[3].legend(frameon=False, fontsize=8)

    fig.suptitle("Ag strip W=800 nm: COMSOL wavelength sweep vs EDP-like estimate")
    fig.tight_layout()
    fig.savefig(results_dir / "comsol_vs_analytic_lambda_sweep_w800.png", dpi=180)
    plt.close(fig)

    print("\n".join(lines))


if __name__ == "__main__":
    main()
