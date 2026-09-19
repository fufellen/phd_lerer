"""Full-vector phase-shifter characteristic and a finite number of PCM levels in the 1D OPA.

Continues simulate_pcm_opa.py with the same array model and metrics.

1. Full-vector (FEM, COMSOL mode analysis) characteristic of the LR-DLSPP section, taken from the
   APEDE-2026 summary (check of 18.07.2026 and recomputation of 18.09.2026): L_pi and attenuation of
   the amorphous and crystalline states. The FEM model has 110 nm of PCM under the gold, so the EIM
   lengths are taken at the same geometry, for the second EIM stage solved in TM polarization (as in
   the article and in VIBR_T) and in TE polarization (the polarization of the vertical field of the
   mode). A device sized by EIM is evaluated with the FEM response; a phase-calibrated device has
   the exact phase and FEM losses.
2. Finite number of PCM levels M: the command phase of a FEM-calibrated device is rounded to a
   multiple of 2 pi / M. Beam metrics are swept over the target angle, and the mean target
   efficiency is checked against the analytical value for a uniformly distributed phase error,
   1/N + (1 - 1/N) sinc^2(pi/M).

Run: python opa_pcm/scripts/fem_and_levels.py   (NumPy and Matplotlib)
"""

from __future__ import annotations

import csv
import math
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import simulate_pcm_opa as S  # noqa: E402


@dataclass(frozen=True)
class FemShifter:
    material: str
    lpi_fem_um: float
    alpha_a_db_um: float
    alpha_c_db_um: float
    lpi_eim_tm_um: float  # EIM, TM on both stages, geometry of the FEM model
    lpi_eim_te_um: float  # EIM, TE on the second stage, same geometry

    @property
    def il_a_db_pi(self) -> float:
        return self.alpha_a_db_um * self.lpi_fem_um

    @property
    def il_c_db_pi(self) -> float:
        return self.alpha_c_db_um * self.lpi_fem_um

    def scale(self, design: str) -> float:
        """True FEM phase / commanded phase for a device sized by the given model."""
        return {"fem": 1.0,
                "eim_tm": self.lpi_eim_tm_um / self.lpi_fem_um,
                "eim_te": self.lpi_eim_te_um / self.lpi_fem_um}[design]


FEM = {
    "GSST": FemShifter("GSST", 0.957, 0.0284, 7.173, 1.045, 0.965),
    "Sb2S3": FemShifter("Sb2S3", 2.824, 0.0196, 0.0208, 4.977, 2.701),
    "Sb2Se3": FemShifter("Sb2Se3", 2.206, 0.0204, 0.0835, 3.042, 2.151),
}
DESIGNS = ("fem", "eim_te", "eim_tm")
LEVELS = (2, 4, 8, 16, 32, 64)


def quantize(q: np.ndarray, levels: int | None) -> np.ndarray:
    if levels is None:
        return q
    step = 2.0 * math.pi / levels
    return (np.round(q / step) % levels) * step


def evaluate_fem(n: int, target_deg: float, sh: FemShifter, design: str,
                 levels: int | None = None, theta_deg=S.THETA_DEG, theta_rad=S.THETA_RAD):
    """Commanded phase q (optionally quantized); the device, sized by `design`, realizes
    scale*q and loses FEM attenuation over its length 2*L_pi(design)."""
    q = quantize(S.target_phases(n, target_deg), levels)
    k = sh.scale(design)
    amp, loss = S.loss_amplitudes(q, sh.il_a_db_pi, sh.il_c_db_pi, length_scale=k)
    power = S.array_power(k * q, amp, theta_rad)
    return S.pattern_metrics(power, amp, target_deg, n, theta_deg, theta_rad), loss, power


def validate() -> None:
    # FEM-calibrated, unquantized: no pointing error at any angle, efficiency = transmission
    for sh in FEM.values():
        for target in (-45.0, 0.0, 17.0, 30.0, 52.0):
            m, _, _ = evaluate_fem(32, target, sh, "fem")
            assert abs(m.pointing_error_deg) <= 0.011, (sh.material, target, m.pointing_error_deg)
    # quantization with lossless elements: the mean over angles follows the uniform-error estimate
    lossless = FemShifter("lossless", 1.0, 0.0, 0.0, 1.0, 1.0)
    angles = np.arange(-60.0, 60.0 + 1e-9, 1.0)
    for m_levels in (4, 8, 16):
        eta = [evaluate_fem(32, a, lossless, "fem", m_levels)[0].eta_target for a in angles]
        x = math.pi / m_levels
        analytic = 1.0 / 32 + (1.0 - 1.0 / 32) * (math.sin(x) / x) ** 2
        assert abs(float(np.mean(eta)) - analytic) < 0.02, (m_levels, np.mean(eta), analytic)
    print("[validate] FEM-calibrated beam points exactly; quantization loss follows sinc^2(pi/M)")


@dataclass
class Row:
    material: str
    N: int
    target_deg: float
    design: str
    levels: int
    phase_scale: float
    mean_loss_db: float
    peak_angle_deg: float
    pointing_error_deg: float
    beamwidth_3db_deg: float
    sll_db: float
    eta_tx: float
    eta_target: float
    eta_coherent: float


def row(sh, n, target, design, levels, m, loss) -> Row:
    return Row(sh.material, n, target, design, 0 if levels is None else levels, sh.scale(design),
               float(np.mean(loss)), m.peak_angle_deg, m.pointing_error_deg, m.beamwidth_3db_deg,
               m.sll_db, m.eta_tx, m.eta_target, m.eta_coherent)


def write(rows: list[Row], name: str) -> None:
    path = S.RESULTS / name
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(asdict(rows[0]).keys()))
        w.writeheader()
        w.writerows(asdict(r) for r in rows)
    print(f"[out] {path}")


def fem_metrics() -> list[Row]:
    rows = []
    for n in S.N_VALUES:
        for sh in FEM.values():
            for design in DESIGNS:
                m, loss, _ = evaluate_fem(n, S.TARGET_DEG, sh, design)
                rows.append(row(sh, n, S.TARGET_DEG, design, None, m, loss))
    write(rows, "fem_metrics.csv")
    return rows


def fem_sweep() -> list[Row]:
    rows = []
    for target in np.arange(-60.0, 60.0 + 1e-9, 5.0):
        for sh in FEM.values():
            for design in DESIGNS:
                m, loss, _ = evaluate_fem(32, float(target), sh, design)
                rows.append(row(sh, 32, float(target), design, None, m, loss))
    write(rows, "fem_steering_sweep_n32.csv")
    return rows


def levels_sweep() -> list[dict]:
    """Per material and M: statistics over target angles -60..60 deg (1 deg) for N = 32."""
    angles = np.arange(-60.0, 60.0 + 1e-9, 1.0)
    theta_deg = np.linspace(-90.0, 90.0, 7201)
    theta_rad = np.deg2rad(theta_deg)
    out = []
    for sh in FEM.values():
        base = [evaluate_fem(32, a, sh, "fem", None, theta_deg, theta_rad)[0] for a in angles]
        for m_levels in LEVELS:
            ms = [evaluate_fem(32, a, sh, "fem", m_levels, theta_deg, theta_rad)[0] for a in angles]
            ratio = np.array([q.eta_target / b.eta_target for q, b in zip(ms, base)])
            x = math.pi / m_levels
            out.append({
                "material": sh.material, "N": 32, "levels": m_levels,
                "eta_ratio_mean": float(np.mean(ratio)),
                "eta_ratio_min": float(np.min(ratio)),
                "eta_ratio_analytic": 1.0 / 32 + (1.0 - 1.0 / 32) * (math.sin(x) / x) ** 2,
                "sll_worst_db": float(max(q.sll_db for q in ms)),
                "sll_unquantized_worst_db": float(max(b.sll_db for b in base)),
                "pointing_error_max_deg": float(max(abs(q.pointing_error_deg) for q in ms)),
            })
    path = S.RESULTS / "pcm_levels_n32.csv"
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)
    print(f"[out] {path}")
    return out


def plot(sweep: list[Row], lv: list[dict]) -> None:
    """(a), (b): device sized by EIM with TM or TE on the second stage, FEM response, N = 32;
    (c), (d): phase quantization of the FEM-calibrated low-loss shifters (GSST is left out: its
    loss law changes with the crystalline fraction and masks the quantization)."""
    fig, axes = plt.subplots(2, 2, figsize=(10.4, 7.2))
    styles = {"GSST": "-", "Sb2S3": "--", "Sb2Se3": ":"}
    names = {"GSST": "GSST", "Sb2S3": "Sb$_2$S$_3$", "Sb2Se3": "Sb$_2$Se$_3$"}
    for design, color, lw, label in (("eim_tm", "0.1", 1.4, "TM"), ("eim_te", "0.6", 2.2, "TE")):
        for mat in FEM:
            sel = [r for r in sweep if r.material == mat and r.design == design]
            x = [r.target_deg for r in sel]
            axes[0, 0].plot(x, [r.eta_target for r in sel], styles[mat], color=color, lw=lw,
                            label=f"{names[mat]}, второй этап ЭДП в {label}")
            axes[0, 1].plot(x, [r.sll_db for r in sel], styles[mat], color=color, lw=lw)
    axes[0, 0].set_ylabel("Норм. мощность в целевом угле")
    axes[0, 0].set_title("(а) проект по ЭДП, ответ по МКЭ, N = 32", fontsize=10)
    axes[0, 0].legend(fontsize=7, loc="center left")
    axes[0, 1].set_ylabel("Макс. боковой лепесток, дБ")
    axes[0, 1].set_title("(б) то же, наибольший боковой лепесток", fontsize=10)
    for ax in axes[0]:
        ax.set_xlabel("Целевой угол, град.")
    low = [m for m in FEM if m != "GSST"]
    for mat in low:
        sel = [d for d in lv if d["material"] == mat]
        axes[1, 0].plot([d["levels"] for d in sel],
                        [10 * math.log10(d["eta_ratio_mean"]) for d in sel],
                        styles[mat], color="0.1", marker="o", ms=3,
                        label=f"{names[mat]}, среднее по углам")
        axes[1, 1].plot([d["levels"] for d in sel], [d["sll_worst_db"] for d in sel],
                        styles[mat], color="0.1", marker="o", ms=3, label=names[mat])
    ms = np.array(LEVELS, dtype=float)
    analytic = 1.0 / 32 + (1.0 - 1.0 / 32) * (np.sin(np.pi / ms) / (np.pi / ms)) ** 2
    axes[1, 0].plot(ms, 10 * np.log10(analytic), color="0.6", lw=4, alpha=0.6, zorder=0,
                    label="оценка $1/N + (1-1/N)\\,\\mathrm{sinc}^2(\\pi/M)$")
    axes[1, 0].set_ylabel("Изменение мощности в целевом угле, дБ")
    axes[1, 0].set_title("(в) квантование фазы, N = 32, углы от −60° до 60°", fontsize=10)
    axes[1, 0].legend(fontsize=7)
    axes[1, 1].axhline(-13.26, color="0.6", lw=1.0, ls="--",
                       label="равномерная решётка без квантования")
    axes[1, 1].set_ylabel("Наибольший боковой лепесток по углам, дБ")
    axes[1, 1].set_title("(г) худший боковой лепесток", fontsize=10)
    axes[1, 1].legend(fontsize=7)
    for ax in axes[1]:
        ax.set_xscale("log", base=2)
        ax.set_xlabel("Число уровней PCM, $M$")
    for ax in axes.flat:
        ax.grid(alpha=0.25)
    fig.tight_layout()
    path = S.RESULTS / "fem_levels_n32.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f"[out] {path}")


def summary(metrics: list[Row], lv: list[dict]) -> None:
    print("\n=== N = 32, target 30 deg, FEM response ===")
    for r in metrics:
        if r.N == 32:
            print(f"{r.material:7s} {r.design:7s} scale={r.phase_scale:.3f} "
                  f"peak={r.peak_angle_deg:+7.2f} err={r.pointing_error_deg:+7.2f} "
                  f"SLL={r.sll_db:7.2f} eta_tx={r.eta_tx:.3f} eta_target={r.eta_target:.3f}")
    print("\n=== PCM levels, N = 32, FEM-calibrated, statistics over -60..60 deg ===")
    for d in lv:
        print(f"{d['material']:7s} M={d['levels']:3d} eta_ratio mean={d['eta_ratio_mean']:.4f} "
              f"(analytic {d['eta_ratio_analytic']:.4f}) min={d['eta_ratio_min']:.4f} "
              f"SLL worst={d['sll_worst_db']:.2f} dB (unquantized {d['sll_unquantized_worst_db']:.2f}) "
              f"max |pointing err|={d['pointing_error_max_deg']:.2f}")


def main() -> None:
    validate()
    metrics = fem_metrics()
    sweep = fem_sweep()
    lv = levels_sweep()
    plot(sweep, lv)
    summary(metrics, lv)


if __name__ == "__main__":
    main()
