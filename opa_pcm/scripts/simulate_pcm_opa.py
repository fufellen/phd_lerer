"""Propagate PCM phase-shifter errors to a 1D optical phased array.

The array-factor model uses isotropic elements and half-wavelength pitch to
isolate phase and amplitude errors. The strict planar phase-shifter model is a
control limit, not a full-vector truth.
"""

from __future__ import annotations

import csv
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

LAMBDA_UM = 1.550
D_OVER_LAMBDA = 0.5
KD = 2.0 * math.pi * D_OVER_LAMBDA
TARGET_DEG = 30.0
N_VALUES = (16, 32, 64)
THETA_DEG = np.linspace(-90.0, 90.0, 18001)
THETA_RAD = np.deg2rad(THETA_DEG)


@dataclass(frozen=True)
class PhaseShifter:
    material: str
    lpi_eim_um: float
    lpi_planar_um: float
    il_a_eim_db_pi: float
    il_c_eim_db_pi: float
    il_a_planar_db_pi: float
    il_c_planar_db_pi: float
    status: str = "ok"

    @property
    def phase_scale(self) -> float:
        """Planar phase / designed phase for a device sized by EIM."""
        return self.lpi_eim_um / self.lpi_planar_um


PHASE_SHIFTERS = {
    "GSST": PhaseShifter(
        "GSST", 1.0316, 1.0532, 0.0370, 7.5759, 0.1864, 4.4566
    ),
    "Sb2S3": PhaseShifter(
        "Sb2S3",
        4.8132,
        1.7859,
        0.4259,
        0.2340,
        1.0819,
        0.4767,
        status="EIM buffer slice below cutoff",
    ),
    "Sb2Se3": PhaseShifter(
        "Sb2Se3", 3.0082, 1.8208, 0.1529, 0.0410, 0.5131, 0.1003
    ),
}


@dataclass
class PatternMetrics:
    peak_angle_deg: float
    pointing_error_deg: float
    beamwidth_3db_deg: float
    sll_db: float
    eta_tx: float
    eta_peak: float
    eta_target: float
    eta_coherent: float
    main_lobe_fraction: float


def target_phases(n_elements: int, target_deg: float) -> np.ndarray:
    n = np.arange(n_elements, dtype=float)
    unwrapped = -n * KD * math.sin(math.radians(target_deg))
    return np.mod(unwrapped, 2.0 * math.pi)


def array_power(
    phases: np.ndarray,
    amplitudes: np.ndarray,
    theta_rad: np.ndarray = THETA_RAD,
) -> np.ndarray:
    n = np.arange(phases.size, dtype=float)
    steering = np.exp(
        1j * (np.outer(np.sin(theta_rad), n * KD) + phases[np.newaxis, :])
    )
    field = steering @ amplitudes
    return np.abs(field) ** 2


def loss_amplitudes(
    wrapped_phases: np.ndarray,
    il_a_db_pi: float,
    il_c_db_pi: float,
    length_scale: float = 1.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Loss of a segmented 0..2pi shifter with total length 2*L_pi."""
    fraction_c = wrapped_phases / (2.0 * math.pi)
    loss_db = 2.0 * length_scale * (
        (1.0 - fraction_c) * il_a_db_pi + fraction_c * il_c_db_pi
    )
    return 10.0 ** (-loss_db / 20.0), loss_db


def scenario_inputs(
    n_elements: int,
    target_deg: float,
    shifter: PhaseShifter | None,
    scenario: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    q = target_phases(n_elements, target_deg)
    if scenario == "ideal":
        return q, np.ones(n_elements), np.zeros(n_elements)
    if shifter is None:
        raise ValueError("A phase shifter is required")
    if scenario == "eim_prediction":
        amp, loss = loss_amplitudes(q, shifter.il_a_eim_db_pi, shifter.il_c_eim_db_pi)
        return q, amp, loss
    if scenario == "planar_calibrated":
        amp, loss = loss_amplitudes(
            q, shifter.il_a_planar_db_pi, shifter.il_c_planar_db_pi
        )
        return q, amp, loss
    if scenario == "phase_scale_only":
        return shifter.phase_scale * q, np.ones(n_elements), np.zeros(n_elements)
    if scenario == "planar_response_to_eim_design":
        amp, loss = loss_amplitudes(
            q,
            shifter.il_a_planar_db_pi,
            shifter.il_c_planar_db_pi,
            length_scale=shifter.phase_scale,
        )
        return shifter.phase_scale * q, amp, loss
    raise ValueError(f"Unknown scenario: {scenario}")


def _beamwidth_3db(theta_deg: np.ndarray, normalized_power: np.ndarray) -> float:
    peak = int(np.argmax(normalized_power))
    above = normalized_power >= 0.5
    left = peak
    right = peak
    while left > 0 and above[left - 1]:
        left -= 1
    while right < above.size - 1 and above[right + 1]:
        right += 1
    return float(theta_deg[right] - theta_deg[left])


def pattern_metrics(
    power: np.ndarray,
    amplitudes: np.ndarray,
    target_deg: float,
    n_elements: int,
    theta_deg: np.ndarray = THETA_DEG,
    theta_rad: np.ndarray = THETA_RAD,
) -> PatternMetrics:
    peak_idx = int(np.argmax(power))
    peak_power = float(power[peak_idx])
    normalized = power / max(peak_power, 1e-300)
    peak_angle = float(theta_deg[peak_idx])
    target_power = float(np.interp(target_deg, theta_deg, power))

    u0 = math.sin(math.radians(target_deg))
    u = np.sin(theta_rad)
    main_mask = np.abs(u - u0) <= 2.0 / n_elements
    sll = 10.0 * math.log10(
        max(float(np.max(normalized[~main_mask])), 1e-300)
    )

    total_integral = float(np.trapz(power, theta_rad))
    main_integral = float(np.trapz(np.where(main_mask, power, 0.0), theta_rad))
    tx = float(np.sum(amplitudes**2) / n_elements)
    denom_coherent = n_elements * float(np.sum(amplitudes**2))

    return PatternMetrics(
        peak_angle_deg=peak_angle,
        pointing_error_deg=peak_angle - target_deg,
        beamwidth_3db_deg=_beamwidth_3db(theta_deg, normalized),
        sll_db=sll,
        eta_tx=tx,
        eta_peak=peak_power / (n_elements**2),
        eta_target=target_power / (n_elements**2),
        eta_coherent=peak_power / max(denom_coherent, 1e-300),
        main_lobe_fraction=main_integral / max(total_integral, 1e-300),
    )


def evaluate(
    n_elements: int,
    target_deg: float,
    shifter: PhaseShifter | None,
    scenario: str,
    theta_deg: np.ndarray = THETA_DEG,
    theta_rad: np.ndarray = THETA_RAD,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, PatternMetrics]:
    phases, amplitudes, losses = scenario_inputs(
        n_elements, target_deg, shifter, scenario
    )
    power = array_power(phases, amplitudes, theta_rad)
    metrics = pattern_metrics(
        power, amplitudes, target_deg, n_elements, theta_deg, theta_rad
    )
    return power, amplitudes, losses, metrics


def validate_model() -> None:
    for n_elements in N_VALUES:
        q = target_phases(n_elements, 0.0)
        power = array_power(q, np.ones(n_elements))
        m = pattern_metrics(power, np.ones(n_elements), 0.0, n_elements)
        expected_bw = math.degrees(2.0 * 0.886 / n_elements)
        assert abs(m.peak_angle_deg) <= 0.011
        assert -13.5 < m.sll_db < -12.8, (n_elements, m.sll_db)
        assert abs(m.beamwidth_3db_deg - expected_bw) / expected_bw < 0.05
        assert abs(m.eta_peak - 1.0) < 1e-12
        assert abs(m.eta_coherent - 1.0) < 1e-12

    n_elements = 32
    n = np.arange(n_elements, dtype=float)
    unwrapped = -n * KD * math.sin(math.radians(TARGET_DEG))
    wrapped = target_phases(n_elements, TARGET_DEG)
    p1 = array_power(unwrapped, np.ones(n_elements))
    p2 = array_power(wrapped, np.ones(n_elements))
    p3 = array_power(wrapped + 0.731, np.ones(n_elements))
    assert np.max(np.abs(p1 - p2)) < 1e-9
    assert np.max(np.abs(p2 - p3)) < 1e-9
    print("[validate] ideal peak, HPBW, SLL and efficiencies OK")
    print("[validate] wrapped phase and common-phase invariance OK")


def write_deterministic_metrics() -> list[dict]:
    rows: list[dict] = []
    scenarios = (
        "ideal",
        "eim_prediction",
        "planar_calibrated",
        "phase_scale_only",
        "planar_response_to_eim_design",
    )
    for n_elements in N_VALUES:
        for material, shifter in PHASE_SHIFTERS.items():
            for scenario in scenarios:
                obj = None if scenario == "ideal" else shifter
                _, _, losses, metrics = evaluate(
                    n_elements, TARGET_DEG, obj, scenario
                )
                rows.append(
                    {
                        "material": material,
                        "N": n_elements,
                        "target_deg": TARGET_DEG,
                        "scenario": scenario,
                        "phase_scale": 1.0
                        if scenario
                        in ("ideal", "eim_prediction", "planar_calibrated")
                        else shifter.phase_scale,
                        "mean_loss_db": float(np.mean(losses)),
                        "min_loss_db": float(np.min(losses)),
                        "max_loss_db": float(np.max(losses)),
                        **asdict(metrics),
                        "status": shifter.status,
                    }
                )
    path = RESULTS / "deterministic_metrics.csv"
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"[out] {path}")
    return rows


def write_steering_sweep() -> list[dict]:
    rows: list[dict] = []
    for target in np.arange(-60.0, 60.0 + 1e-9, 5.0):
        for material, shifter in PHASE_SHIFTERS.items():
            for scenario in (
                "eim_prediction",
                "planar_calibrated",
                "planar_response_to_eim_design",
            ):
                _, _, losses, metrics = evaluate(32, float(target), shifter, scenario)
                rows.append(
                    {
                        "material": material,
                        "N": 32,
                        "target_deg": float(target),
                        "scenario": scenario,
                        "phase_scale": shifter.phase_scale,
                        "mean_loss_db": float(np.mean(losses)),
                        **asdict(metrics),
                        "status": shifter.status,
                    }
                )
    path = RESULTS / "steering_sweep_n32.csv"
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"[out] {path}")
    return rows


def random_phase_tolerance(n_trials: int = 100) -> list[dict]:
    rng = np.random.default_rng(20260710)
    theta_deg = np.linspace(-90.0, 90.0, 1801)
    theta_rad = np.deg2rad(theta_deg)
    rows: list[dict] = []
    for n_elements in N_VALUES:
        q = target_phases(n_elements, TARGET_DEG)
        amplitudes = np.ones(n_elements)
        for sigma_deg in (0.0, 2.0, 5.0, 10.0, 20.0, 30.0):
            values = []
            trials = 1 if sigma_deg == 0.0 else n_trials
            for _ in range(trials):
                noise = rng.normal(0.0, math.radians(sigma_deg), n_elements)
                power = array_power(q + noise, amplitudes, theta_rad)
                values.append(
                    pattern_metrics(
                        power,
                        amplitudes,
                        TARGET_DEG,
                        n_elements,
                        theta_deg,
                        theta_rad,
                    )
                )
            sigma_rad = math.radians(sigma_deg)
            analytic_eta = (
                1.0 / n_elements
                + (1.0 - 1.0 / n_elements) * math.exp(-(sigma_rad**2))
            )
            for metric_name in (
                "pointing_error_deg",
                "sll_db",
                "eta_target",
                "eta_peak",
                "eta_coherent",
            ):
                arr = np.array([getattr(v, metric_name) for v in values], dtype=float)
                rows.append(
                    {
                        "N": n_elements,
                        "target_deg": TARGET_DEG,
                        "sigma_deg": sigma_deg,
                        "metric": metric_name,
                        "mean": float(np.mean(arr)),
                        "std": float(np.std(arr)),
                        "p05": float(np.quantile(arr, 0.05)),
                        "p95": float(np.quantile(arr, 0.95)),
                        "analytic_eta_target": analytic_eta
                        if metric_name == "eta_target"
                        else "",
                    }
                )
            if sigma_deg > 0.0:
                eta_mean = np.mean([v.eta_target for v in values])
                assert abs(eta_mean - analytic_eta) < 0.035
    path = RESULTS / "random_phase_tolerance.csv"
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"[out] {path}")
    print("[validate] Monte Carlo target efficiency agrees with analytic expectation")
    return rows


def plot_patterns() -> None:
    fig, axes = plt.subplots(3, 1, figsize=(8.0, 9.0), sharex=True)
    ideal, _, _, _ = evaluate(32, TARGET_DEG, None, "ideal")
    ideal_db = 10.0 * np.log10(np.maximum(ideal / ideal.max(), 1e-8))
    for ax, (material, shifter) in zip(axes, PHASE_SHIFTERS.items()):
        ax.plot(THETA_DEG, ideal_db, color="black", lw=1.4, label="идеальная")
        for scenario, label, style in (
            ("eim_prediction", "прогноз ЭДП", "--"),
            ("planar_response_to_eim_design", "планарный ответ ЭДП-проекта", "-"),
        ):
            power, _, _, _ = evaluate(32, TARGET_DEG, shifter, scenario)
            db = 10.0 * np.log10(np.maximum(power / power.max(), 1e-8))
            ax.plot(THETA_DEG, db, style, lw=1.3, label=label)
        ax.axvline(TARGET_DEG, color="0.55", lw=0.8, ls=":")
        ax.set_xlim(-90, 90)
        ax.set_ylim(-50, 1)
        ax.set_ylabel("Норм. мощность, дБ")
        ax.set_title(f"{material}, N=32, целевой угол 30°")
        ax.grid(alpha=0.25)
        ax.legend(loc="lower left", fontsize=8)
    axes[-1].set_xlabel("Угол, град.")
    fig.tight_layout()
    path = RESULTS / "patterns_n32_theta30.png"
    fig.savefig(path, dpi=300)
    plt.close(fig)
    print(f"[out] {path}")


def plot_deterministic_metrics(rows: list[dict]) -> None:
    selected = [
        r
        for r in rows
        if r["N"] == 32 and r["scenario"] == "planar_response_to_eim_design"
    ]
    labels = [r["material"] for r in selected]
    x = np.arange(len(labels))
    fig, axes = plt.subplots(2, 2, figsize=(8.4, 6.6))
    data = (
        ("pointing_error_deg", "Ошибка направления, град."),
        ("sll_db", "Макс. боковой лепесток, дБ"),
        ("eta_tx", "Доля прошедшей мощности"),
        ("eta_target", "Норм. мощность в целевом угле"),
    )
    for ax, (key, ylabel) in zip(axes.flat, data):
        vals = [float(r[key]) for r in selected]
        ax.bar(x, vals, color=("#4C78A8", "#F58518", "#54A24B"))
        ax.set_xticks(x, labels)
        ax.set_ylabel(ylabel)
        ax.grid(axis="y", alpha=0.25)
    fig.suptitle("ЭДП-проект в строгом планарном пределе; N=32, 30°")
    fig.tight_layout()
    path = RESULTS / "metrics_n32_theta30.png"
    fig.savefig(path, dpi=300)
    plt.close(fig)
    print(f"[out] {path}")


def plot_random_tolerance(rows: list[dict]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.5))
    for n_elements in N_VALUES:
        eta = [
            r
            for r in rows
            if r["N"] == n_elements and r["metric"] == "eta_target"
        ]
        sll = [
            r for r in rows if r["N"] == n_elements and r["metric"] == "sll_db"
        ]
        axes[0].plot(
            [r["sigma_deg"] for r in eta],
            [r["mean"] for r in eta],
            marker="o",
            label=f"N={n_elements}",
        )
        axes[1].plot(
            [r["sigma_deg"] for r in sll],
            [r["mean"] for r in sll],
            marker="o",
            label=f"N={n_elements}",
        )
    axes[0].set_xlabel("СКО случайной фазы, град.")
    axes[0].set_ylabel("Средняя мощность в целевом угле")
    axes[1].set_xlabel("СКО случайной фазы, град.")
    axes[1].set_ylabel("Средний макс. боковой лепесток, дБ")
    for ax in axes:
        ax.grid(alpha=0.25)
        ax.legend()
    fig.tight_layout()
    path = RESULTS / "random_phase_tolerance.png"
    fig.savefig(path, dpi=300)
    plt.close(fig)
    print(f"[out] {path}")


def plot_steering_sweep(rows: list[dict]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.6))
    for material in PHASE_SHIFTERS:
        selected = [
            r
            for r in rows
            if r["material"] == material
            and r["scenario"] == "planar_response_to_eim_design"
        ]
        axes[0].plot(
            [r["target_deg"] for r in selected],
            [r["pointing_error_deg"] for r in selected],
            marker=".",
            label=material,
        )
        axes[1].plot(
            [r["target_deg"] for r in selected],
            [r["eta_target"] for r in selected],
            marker=".",
            label=material,
        )
    axes[0].axhline(0.0, color="black", lw=0.8)
    axes[0].set_xlabel("Целевой угол, град.")
    axes[0].set_ylabel("Ошибка направления, град.")
    axes[1].set_xlabel("Целевой угол, град.")
    axes[1].set_ylabel("Норм. мощность в целевом угле")
    for ax in axes:
        ax.grid(alpha=0.25)
        ax.legend()
    fig.tight_layout()
    path = RESULTS / "steering_sweep_n32.png"
    fig.savefig(path, dpi=300)
    plt.close(fig)
    print(f"[out] {path}")


def print_primary_summary(rows: list[dict]) -> None:
    print("\n=== N=32, target 30 deg: planar response to EIM design ===")
    for r in rows:
        if r["N"] == 32 and r["scenario"] == "planar_response_to_eim_design":
            print(
                f"{r['material']:7s} phase_scale={r['phase_scale']:.3f} "
                f"peak={r['peak_angle_deg']:+.2f} deg "
                f"err={r['pointing_error_deg']:+.2f} deg "
                f"SLL={r['sll_db']:.2f} dB "
                f"eta_tx={r['eta_tx']:.3f} eta_target={r['eta_target']:.3f}"
            )


def main() -> None:
    validate_model()
    deterministic = write_deterministic_metrics()
    steering = write_steering_sweep()
    random_rows = random_phase_tolerance()
    plot_patterns()
    plot_deterministic_metrics(deterministic)
    plot_random_tolerance(random_rows)
    plot_steering_sweep(steering)
    print_primary_summary(deterministic)


if __name__ == "__main__":
    main()
