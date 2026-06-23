from __future__ import annotations

import csv
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

try:
    from scipy.optimize import root as scipy_root
except ModuleNotFoundError:
    scipy_root = None


BASE_DIR = Path(__file__).resolve().parents[2]
AG_C = BASE_DIR / "Ag_.c"

N_SUB = 1.45
EPS_AIR = 1.0 + 0.0j
EPS_SUB = N_SUB**2 + 0.0j
T_AG_UM = 0.020
W_STRIP_UM = 0.800
LAMBDAS_NM = np.arange(450.0, 801.0, 25.0)
TABLE_LAMBDAS_NM = np.arange(450.0, 801.0, 50.0)


def load_nk_from_c(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    text = path.read_text(encoding="utf-8", errors="ignore")
    pattern = re.compile(
        r"Lam\[(\d+)\]\s*=\s*([0-9.]+);\s*"
        r"n\[\1\]\s*=\s*([0-9.]+);\s*"
        r"k\[\1\]\s*=\s*([0-9.]+);"
    )
    rows = []
    for match in pattern.finditer(text):
        rows.append((float(match.group(2)), float(match.group(3)), float(match.group(4))))
    if not rows:
        raise RuntimeError(f"No optical constants found in {path}")
    rows.sort()
    lam, n, k = (np.array(col, dtype=float) for col in zip(*rows))
    return lam, n, k


AG_LAM_NM, AG_N, AG_K = load_nk_from_c(AG_C)


def ag_nk(lambda_nm: float) -> tuple[float, float]:
    n = float(np.interp(lambda_nm, AG_LAM_NM, AG_N))
    k = float(np.interp(lambda_nm, AG_LAM_NM, AG_K))
    return n, k


def ag_eps(lambda_nm: float, loss_sign: int = +1) -> complex:
    n, k = ag_nk(lambda_nm)
    return complex(n, loss_sign * k) ** 2


def choose_decay_sqrt(z: complex) -> complex:
    value = np.sqrt(z + 0j)
    if value.real < 0:
        value = -value
    if abs(value.real) < 1e-12 and value.imag < 0:
        value = -value
    return value


def planar_tm_residual(neff: complex, lambda_nm: float) -> complex:
    lambda_um = lambda_nm / 1000.0
    k0 = 2.0 * math.pi / lambda_um
    eps = [EPS_AIR, ag_eps(lambda_nm), EPS_SUB]
    gamma = [k0 * choose_decay_sqrt(neff * neff - e) for e in eps]
    q = [gamma[i] / eps[i] for i in range(3)]
    g2 = gamma[1]
    return (q[0] + q[2]) * q[1] * np.cosh(g2 * T_AG_UM) + (
        q[0] * q[2] + q[1] * q[1]
    ) * np.sinh(g2 * T_AG_UM)


def solve_complex(func, guess: complex, args: tuple) -> complex:
    if scipy_root is None:
        return solve_complex_newton(func, guess, args)

    def wrapped(v: np.ndarray) -> np.ndarray:
        value = func(complex(v[0], v[1]), *args)
        return np.array([value.real, value.imag])

    sol = scipy_root(wrapped, np.array([guess.real, guess.imag]), method="hybr", tol=1e-11)
    if not sol.success:
        raise RuntimeError(f"Root solve failed near {guess}: {sol.message}")
    return complex(sol.x[0], sol.x[1])


def solve_complex_newton(func, guess: complex, args: tuple) -> complex:
    z = complex(guess)
    for _ in range(80):
        value = func(z, *args)
        f0 = np.array([value.real, value.imag], dtype=float)
        norm0 = float(np.linalg.norm(f0, ord=2))
        if norm0 < 1e-11:
            return z

        hx = max(1e-7, abs(z.real) * 1e-7)
        hy = max(1e-7, abs(z.imag) * 1e-7)
        fx = func(z + hx, *args)
        fy = func(z + 1j * hy, *args)
        jac = np.array(
            [
                [(fx.real - value.real) / hx, (fy.real - value.real) / hy],
                [(fx.imag - value.imag) / hx, (fy.imag - value.imag) / hy],
            ],
            dtype=float,
        )

        try:
            step = np.linalg.solve(jac, -f0)
        except np.linalg.LinAlgError:
            step = np.linalg.lstsq(jac, -f0, rcond=None)[0]

        accepted = False
        scale = 1.0
        for _ in range(14):
            candidate = z + scale * complex(step[0], step[1])
            candidate_value = func(candidate, *args)
            candidate_norm = abs(candidate_value)
            if candidate_norm < norm0:
                z = candidate
                accepted = True
                break
            scale *= 0.5

        if not accepted:
            z += 0.1 * complex(step[0], step[1])

        if abs(scale * complex(step[0], step[1])) < 1e-12 * max(1.0, abs(z)):
            return z

    raise RuntimeError(f"Newton root solve failed near {guess}: residual={func(z, *args)}")


def try_solve_complex(func, guess: complex, args: tuple) -> complex | None:
    try:
        return solve_complex(func, guess, args)
    except RuntimeError:
        return None


def solve_planar_branch(lambdas_nm: np.ndarray) -> dict[float, complex]:
    out: dict[float, complex] = {}
    guess = 2.40 + 0.20j
    for lam in lambdas_nm:
        neff = solve_complex(planar_tm_residual, guess, (float(lam),))
        if neff.imag < 0:
            neff = neff.conjugate()
        out[float(lam)] = neff
        guess = neff
    return out


def horizontal_tm_residual(
    neff: complex, lambda_nm: float, eps_core: complex, width_um: float = W_STRIP_UM
) -> complex:
    lambda_um = lambda_nm / 1000.0
    k0 = 2.0 * math.pi / lambda_um
    half_width = width_um / 2.0
    u = k0 * np.sqrt(eps_core - neff * neff + 0j)
    alpha = k0 * choose_decay_sqrt(neff * neff - EPS_SUB)
    return (u / eps_core) * np.tan(u * half_width) - alpha / EPS_SUB


def solve_strip_branch(
    planar: dict[float, complex], width_um: float = W_STRIP_UM
) -> dict[float, complex]:
    out: dict[float, complex] = {}
    previous: complex | None = None
    for lam, n_planar in planar.items():
        eps_core = n_planar * n_planar
        guesses = []
        if previous is not None:
            guesses.extend(
                [
                    previous,
                    previous + 0.02,
                    previous - 0.02,
                    previous + 0.02j,
                    previous - 0.02j,
                ]
            )
        guesses.extend(
            [
                n_planar - 0.0005,
                n_planar - 0.0010,
                n_planar - 0.0020,
                n_planar - 0.005,
                n_planar - 0.010,
                complex(n_planar.real - 0.005, n_planar.imag - 0.02),
                complex(n_planar.real - 0.010, n_planar.imag - 0.02),
                complex(n_planar.real - 0.005, n_planar.imag + 0.02),
                complex(n_planar.real - 0.010, n_planar.imag + 0.02),
                n_planar - 0.02,
                n_planar - 0.05,
                n_planar - 0.10,
                complex(n_planar.real - 0.12, n_planar.imag),
                complex(n_planar.real - 0.20, n_planar.imag),
                complex(max(N_SUB + 0.05, n_planar.real - 0.35), n_planar.imag),
            ]
        )
        roots: list[complex] = []
        for guess in guesses:
            candidate = try_solve_complex(
                horizontal_tm_residual, guess, (float(lam), eps_core, width_um)
            )
            if candidate is None:
                continue
            if candidate.imag < 0:
                candidate = candidate.conjugate()
            if candidate.real <= N_SUB or candidate.real > n_planar.real + 1e-6:
                continue
            if not any(abs(candidate - existing) < 1e-6 for existing in roots):
                roots.append(candidate)
        if not roots:
            raise RuntimeError(f"No guided strip root found at lambda={lam} nm")
        neff = max(roots, key=lambda z: z.real)
        if neff.imag < 0:
            neff = neff.conjugate()
        out[float(lam)] = neff
        previous = neff
    return out


def attenuation_ratio(neff: complex) -> float:
    return abs(neff.imag / neff.real)


def write_csv(path: Path, planar: dict[float, complex], strip: dict[float, complex]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "lambda_nm",
                "ag_n",
                "ag_k",
                "ag_eps_real",
                "ag_eps_imag",
                "neff_infinite_real",
                "neff_infinite_imag",
                "ratio_infinite",
                "neff_w800_real",
                "neff_w800_imag",
                "ratio_w800",
            ]
        )
        for lam in planar:
            n, k = ag_nk(lam)
            eps = ag_eps(lam)
            writer.writerow(
                [
                    f"{lam:.0f}",
                    f"{n:.6g}",
                    f"{k:.6g}",
                    f"{eps.real:.9g}",
                    f"{eps.imag:.9g}",
                    f"{planar[lam].real:.9g}",
                    f"{planar[lam].imag:.9g}",
                    f"{attenuation_ratio(planar[lam]):.9g}",
                    f"{strip[lam].real:.9g}",
                    f"{strip[lam].imag:.9g}",
                    f"{attenuation_ratio(strip[lam]):.9g}",
                ]
            )


def write_summary(path: Path, planar: dict[float, complex], strip: dict[float, complex]) -> None:
    lines = []
    lines.append("Analytic repeat of Lerer's metal strip graph")
    lines.append("Geometry: air | 20 nm Ag | substrate n=1.45; finite strip W=800 nm.")
    lines.append("Silver optical constants are parsed from Ag_.c. COMSOL/RF loss sign eps=(n+i*k)^2 is used.")
    lines.append("")
    lines.append("lambda_nm  n_inf_real  n_inf_imag  ratio_inf  n_w800_real  n_w800_imag  ratio_w800")
    for lam in TABLE_LAMBDAS_NM:
        p = planar[float(lam)]
        s = strip[float(lam)]
        lines.append(
            f"{lam:8.0f}  {p.real:10.6f}  {p.imag:10.6f}  {attenuation_ratio(p):9.6f}  "
            f"{s.real:11.6f}  {s.imag:11.6f}  {attenuation_ratio(s):10.6f}"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def plot(path: Path, planar: dict[float, complex], strip: dict[float, complex]) -> None:
    lam = np.array(list(planar.keys()))
    n_inf = np.array([planar[x] for x in lam])
    n_strip = np.array([strip[x] for x in lam])

    fig, (ax_loss, ax_neff) = plt.subplots(
        2, 1, figsize=(6.2, 6.0), sharex=True, gridspec_kw={"height_ratios": [1, 1.45]}
    )
    ax_loss.plot(lam, [attenuation_ratio(x) for x in n_strip], "s-", ms=3, label="W=800 nm")
    ax_loss.plot(lam, [attenuation_ratio(x) for x in n_inf], "o-", ms=3, label="W=infinity")
    ax_loss.set_ylabel("n''/n'")
    ax_loss.grid(True, color="0.82")
    ax_loss.legend(frameon=False, fontsize=9)

    ax_neff.plot(lam, n_strip.real, "s-", ms=3, label="W=800 nm")
    ax_neff.plot(lam, n_inf.real, "o-", ms=3, label="W=infinity")
    ax_neff.set_ylabel("n'")
    ax_neff.set_xlabel("lambda, nm")
    ax_neff.grid(True, color="0.82")
    ax_neff.legend(frameon=False, fontsize=9)

    fig.suptitle("Ag strip on n=1.45 substrate, t=20 nm")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> None:
    out_dir = Path(__file__).resolve().parents[1] / "results"
    out_dir.mkdir(parents=True, exist_ok=True)
    planar = solve_planar_branch(LAMBDAS_NM)
    strip = solve_strip_branch(planar)
    write_csv(out_dir / "analytic_metal_strip_w800.csv", planar, strip)
    write_summary(out_dir / "analytic_metal_strip_w800_summary.txt", planar, strip)
    plot(out_dir / "analytic_metal_strip_w800.png", planar, strip)


if __name__ == "__main__":
    main()
