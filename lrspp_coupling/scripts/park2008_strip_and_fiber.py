"""Аналитическое воспроизведение Park et al., Opt. Commun. 281, 2057 (2008).

Статья: полосковые LR-SPP-волноводы в симметричном полимере на 1.31 и 1.55 мкм.
DOI 10.1016/j.optcom.2007.10.115.

Что воспроизводится (числа статьи, все измеренные - методом cutback):
  потери на грань, 1.55 мкм, Au 14.5 нм: 4.45 / 3.20 / 1.17 / 0.37 дБ
    при ширинах 1.5 / 2.0 / 3.0 / 5.0 мкм
  потери на грань, 1.31 мкм: плато около 0.6 дБ при ширинах 3.5-4.5 мкм
  потери распространения, 1.55 мкм: 1.40 дБ/см при 2 мкм, 4.30 дБ/см при 5 мкм
  потери распространения, 1.31 мкм: 2.05 дБ/см при 1.5 мкм, 11.15 дБ/см при 5 мкм
  расчёт авторов методом линий, ширина 5 мкм, 1.55 мкм: 0.2 / 1.15 / 3.85 / 6.2
    дБ/см при толщинах 10 / 14 / 18 / 20 нм

Главная новая величина по сравнению с остальными скриптами - потери на грань.
Они считаются как перекрытие двумерного профиля моды полоски с гауссовым полем
одномодового волокна, то есть ровно та физика, ради которой в обзоре собиралась
таблица ввода излучения.

Запуск:
    python lrspp_coupling/scripts/park2008_strip_and_fiber.py
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from slabmodes import Stack, materials, propagation_loss_db_per_cm, solve_mode, solve_strip  # noqa: E402
from slabmodes.coupling import (  # noqa: E402
    coupling_loss_db,
    gaussian_field,
    overlap_2d,
    strip_mode_field,
)

OUT = Path(__file__).resolve().parents[1] / "results"
OUT.mkdir(parents=True, exist_ok=True)

T_AU_UM = 0.0145  # толщина золота в измерениях рис. 4 статьи

# Длина волны -> (eps золота, eps полимера, диаметр модового пятна волокна)
CASES = {
    1.55: (materials.EPS_AU_1550, materials.EPS_ZPU450, 10.4),
    1.31: (materials.EPS_AU_1310, materials.EPS_ZPU450_1310, 9.2),
}

WIDTHS = np.array([1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0])

# Считанные с рис. 4 статьи значения (эксперимент)
PAPER_FACET = {
    1.55: {1.5: 4.45, 2.0: 3.20, 2.5: 1.92, 3.0: 1.17, 3.5: 0.75, 4.0: 0.57, 4.5: 0.50, 5.0: 0.37},
    1.31: {1.5: 1.63, 2.0: 1.10, 2.5: 0.80, 3.0: 0.70, 3.5: 0.58, 4.0: 0.58, 4.5: 0.58, 5.0: 0.62},
}
PAPER_PROP = {
    1.55: {1.5: 1.45, 2.0: 1.40, 2.5: 1.95, 3.0: 2.35, 3.5: 2.85, 4.0: 3.25, 4.5: 3.75, 5.0: 4.30},
    1.31: {1.5: 2.05, 2.0: 2.90, 2.5: 4.40, 3.0: 5.95, 3.5: 7.50, 4.0: 8.95, 4.5: 10.10, 5.0: 11.15},
}
# Расчёт авторов методом линий, ширина 5 мкм, 1.55 мкм: толщина -> дБ/см
PAPER_CALC_THICKNESS = {0.010: 0.2, 0.014: 1.15, 0.018: 3.85, 0.020: 6.2, 0.024: 13.2}


def stack_for(lam: float, thickness_um: float = T_AU_UM) -> Stack:
    eps_au, eps_pol, _ = CASES[lam]
    return Stack(eps=(eps_pol, eps_au, eps_pol), thickness=(thickness_um,))


def facet_loss(lam: float, width_um: float, n_planar: complex, n_strip: complex) -> float:
    """Потери стыковки одномодового волокна с полоской, дБ на грань."""
    _, eps_pol, mfd = CASES[lam]
    k0 = materials.k0_from_lambda(lam)
    x = np.linspace(-22.0, 22.0, 901)
    z = np.linspace(-22.0, 22.0, 901)
    mode = strip_mode_field(n_planar, n_strip, stack_for(lam), k0, width_um, eps_pol, x, z)
    fiber = gaussian_field(x, z, mfd)
    return coupling_loss_db(overlap_2d(mode, fiber, x, z))


def width_sweep(lam: float) -> list[dict]:
    _, eps_pol, _ = CASES[lam]
    k0 = materials.k0_from_lambda(lam)
    n_planar = solve_mode(stack_for(lam), k0, complex(np.sqrt(eps_pol).real + 1.5e-3, 2e-5))
    rows = []
    for w in WIDTHS:
        try:
            n_strip = solve_strip(k0, n_planar**2, eps_pol, float(w))
        except RuntimeError:
            rows.append({"width": float(w), "strip": None})
            continue
        rows.append({
            "width": float(w),
            "planar": n_planar,
            "strip": n_strip,
            "prop_db_cm": propagation_loss_db_per_cm(n_strip, lam),
            "facet_db": facet_loss(lam, float(w), n_planar, n_strip),
        })
    return rows


def thickness_sweep(lam: float = 1.55, width_um: float = 5.0) -> list[dict]:
    _, eps_pol, _ = CASES[lam]
    k0 = materials.k0_from_lambda(lam)
    rows = []
    for t in (0.010, 0.012, 0.014, 0.016, 0.018, 0.020, 0.022, 0.024):
        stack = Stack(eps=(eps_pol, CASES[lam][0], eps_pol), thickness=(t,))
        n_planar = solve_mode(stack, k0, complex(np.sqrt(eps_pol).real + 2e-3, 5e-5))
        try:
            n_strip = solve_strip(k0, n_planar**2, eps_pol, width_um)
        except RuntimeError:
            n_strip = None
        rows.append({
            "thickness": t,
            "planar_db_cm": propagation_loss_db_per_cm(n_planar, lam),
            "strip_db_cm": propagation_loss_db_per_cm(n_strip, lam) if n_strip else None,
        })
    return rows


def make_plots(sweeps: dict[float, list[dict]], thick: list[dict]) -> None:
    fig, axs = plt.subplots(1, 3, figsize=(15.6, 4.9))

    ax = axs[0]
    for lam, color in ((1.55, "#0B6E99"), (1.31, "#B23A48")):
        rows = [r for r in sweeps[lam] if r.get("strip")]
        ax.plot([r["width"] for r in rows], [r["facet_db"] for r in rows], "o-", ms=4,
                color=color, label=f"расчёт, {lam:g} мкм")
        ref = PAPER_FACET[lam]
        ax.plot(list(ref), list(ref.values()), "*", ms=13, color=color,
                label=f"измерено, {lam:g} мкм")
    ax.set_xlabel("ширина полоски, мкм")
    ax.set_ylabel("потери на грань, дБ")
    ax.set_title("Стыковка с одномодовым волокном")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, frameon=False)

    ax = axs[1]
    for lam, color in ((1.55, "#0B6E99"), (1.31, "#B23A48")):
        rows = [r for r in sweeps[lam] if r.get("strip")]
        ax.plot([r["width"] for r in rows], [r["prop_db_cm"] for r in rows], "o-", ms=4,
                color=color, label=f"расчёт, {lam:g} мкм")
        ref = PAPER_PROP[lam]
        ax.plot(list(ref), list(ref.values()), "*", ms=13, color=color,
                label=f"измерено, {lam:g} мкм")
    ax.set_xlabel("ширина полоски, мкм")
    ax.set_ylabel("потери распространения, дБ/см")
    ax.set_title("Потери против ширины")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, frameon=False)

    ax = axs[2]
    t_nm = [r["thickness"] * 1000 for r in thick]
    ax.semilogy(t_nm, [r["planar_db_cm"] for r in thick], "o-", ms=4, label="расчёт, плёнка")
    strip = [(r["thickness"] * 1000, r["strip_db_cm"]) for r in thick if r["strip_db_cm"]]
    if strip:
        ax.semilogy([s[0] for s in strip], [s[1] for s in strip], "s-", ms=4,
                    label="расчёт, полоска 5 мкм")
    ax.semilogy([t * 1000 for t in PAPER_CALC_THICKNESS],
                list(PAPER_CALC_THICKNESS.values()), "k*", ms=13,
                label="расчёт авторов, метод линий")
    ax.set_xlabel("толщина золота, нм")
    ax.set_ylabel("потери, дБ/см")
    ax.set_title("Потери против толщины, 1.55 мкм")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8, frameon=False)

    fig.suptitle("Park 2008: полосковые LR-SPP-волноводы и ввод из волокна", fontsize=12)
    fig.tight_layout()
    fig.savefig(OUT / "park2008_strip.png", dpi=170)
    plt.close(fig)


def main() -> int:
    sweeps = {lam: width_sweep(lam) for lam in CASES}
    thick = thickness_sweep()
    make_plots(sweeps, thick)

    lines: list[str] = []
    add = lines.append
    add("Park et al., Opt. Commun. 281, 2057 (2008) - аналитическое воспроизведение")
    add(f"Полоска золота {T_AU_UM * 1000:.1f} нм в симметричном полимере ZPU450")
    add("")
    for lam in (1.55, 1.31):
        _, eps_pol, mfd = CASES[lam]
        rows = [r for r in sweeps[lam] if r.get("strip")]
        add(f"Длина волны {lam:g} мкм, показатель полимера {np.sqrt(eps_pol).real:.4f}, "
            f"пятно волокна {mfd:g} мкм")
        add(f"  планарный предел n_eff = {rows[0]['planar'].real:.6f}"
            f"   потери {propagation_loss_db_per_cm(rows[0]['planar'], lam):7.3f} дБ/см")
        add("  W, мкм    n_eff    потери расч.  потери изм.   грань расч.  грань изм.")
        for r in rows:
            w = r["width"]
            pr = PAPER_PROP[lam].get(w)
            fr = PAPER_FACET[lam].get(w)
            add(f"  {w:5.1f}   {r['strip'].real:.6f}   {r['prop_db_cm']:10.2f}   "
                f"{(f'{pr:9.2f}' if pr else '        -')}   {r['facet_db']:10.2f}   "
                f"{(f'{fr:8.2f}' if fr else '       -')}")
        add("")

    add("Потери против толщины при ширине 5 мкм, 1.55 мкм")
    add("  t, нм   плёнка расч.  полоска расч.  расчёт авторов")
    for r in thick:
        ref = PAPER_CALC_THICKNESS.get(r["thickness"])
        sv = f"{r['strip_db_cm']:12.2f}" if r["strip_db_cm"] else "           -"
        add(f"  {r['thickness'] * 1000:5.0f}   {r['planar_db_cm']:11.2f}   {sv}   "
            f"{(f'{ref:13.2f}' if ref else '            -')}")

    with (OUT / "park2008_strip.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["lambda_um", "width_um", "neff_real", "neff_imag",
                    "prop_db_cm_calc", "prop_db_cm_meas", "facet_db_calc", "facet_db_meas"])
        for lam in CASES:
            for r in sweeps[lam]:
                if not r.get("strip"):
                    continue
                w.writerow([f"{lam:g}", f"{r['width']:.1f}",
                            f"{r['strip'].real:.9f}", f"{r['strip'].imag:.6e}",
                            f"{r['prop_db_cm']:.4f}",
                            PAPER_PROP[lam].get(r["width"], ""),
                            f"{r['facet_db']:.4f}",
                            PAPER_FACET[lam].get(r["width"], "")])

    (OUT / "park2008_summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    for line in lines:
        print(line)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
