"""Аналитическое воспроизведение Park et al., Opt. Commun. 282, 4513 (2009).

Статья: вертикальный направленный ответвитель между LR-SPP-волноводом (полоска
золота) и обычным диэлектрическим волноводом, разнесёнными по вертикали.
DOI 10.1016/j.optcom.2009.08.038.

Что воспроизводится и с чем сверяется (числа статьи):
  n_eff обоих волноводов               ~1.4307
  длина связи при зазоре 5 мкм         680 мкм (расчёт), 600 мкм (эксперимент)
  длина связи при 6 / 8 / 10 мкм       880 / ~1480 / ~2450 мкм
  отсечка нижней супермоды             зазор ~4.8 мкм
  максимальная перекачка               98 % (расчёт), 86 % (эксперимент)
  затухание устройства                 4.1 дБ/см (расчёт), 7.6 дБ/см (эксперимент)
  потери одиночного LR-SPP             ~9 дБ/см (расчёт), 13 дБ/см (эксперимент)

Метод и главная тонкость. Вертикальная задача решается точно методом матрицы
переноса для пятислойной структуры. Но в планарном пределе волноводы НЕ
синхронизованы: их эффективные показатели различаются на 2.2e-3, что больше
связи, поэтому расщепление супермод определяется расстройкой и длина связи
получается вдвое короче наблюдаемой. Синхронизм в приборе возникает только
после сужения полосок по ширине. Поэтому расчёт разделён на два шага:

  1) из точных планарных супермод извлекается коэффициент связи kappa - он
     задаётся перекрытием затухающих хвостов поперёк зазора и от ширины полосок
     почти не зависит;
  2) расстройка берётся реальная, по эффективным показателям полосок конечной
     ширины (метод эффективного показателя).

Запуск:
    python lrspp_coupling/scripts/park2009_vertical_coupler.py
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

from slabmodes import (  # noqa: E402
    CoupledPair,
    Stack,
    extract_kappa,
    materials,
    mode_fields,
    propagation_loss_db_per_cm,
    solve_mode,
    solve_strip,
)

OUT = Path(__file__).resolve().parents[1] / "results"
OUT.mkdir(parents=True, exist_ok=True)

LAMBDA_UM = 1.55
K0 = materials.k0_from_lambda(LAMBDA_UM)

EPS_CLAD = materials.EPS_ZPU13_430
EPS_CORE = materials.EPS_ZPU13_440
EPS_AU = materials.EPS_AU_1550

T_AU_UM = 0.020
W_AU_UM = 6.0
T_CORE_UM = 2.8
W_CORE_UM = 3.0
GAP_NOMINAL_UM = 5.0

PAPER_LC = {5.0: 680.0, 6.0: 880.0, 8.0: 1480.0, 10.0: 2450.0}
PAPER_NEFF = 1.4307
PAPER_CUTOFF_GAP = 4.8
PAPER_MAX_TRANSFER = 0.98
PAPER_DEVICE_LOSS = 4.1
PAPER_LRSPP_LOSS = 9.0


def film_stack() -> Stack:
    return Stack(eps=(EPS_CLAD, EPS_AU, EPS_CLAD), thickness=(T_AU_UM,),
                 names=("полимер", "Au", "полимер"))


def slab_stack() -> Stack:
    return Stack(eps=(EPS_CLAD, EPS_CORE, EPS_CLAD), thickness=(T_CORE_UM,),
                 names=("полимер", "ZPU13-440", "полимер"))


def coupled_stack(gap_um: float) -> Stack:
    return Stack(
        eps=(EPS_CLAD, EPS_AU, EPS_CLAD, EPS_CORE, EPS_CLAD),
        thickness=(T_AU_UM, gap_um, T_CORE_UM),
        names=("полимер", "Au", "зазор", "ZPU13-440", "полимер"),
    )


def isolated_modes() -> dict[str, complex]:
    """Планарные и полосковые эффективные показатели обоих волноводов."""
    n_film = solve_mode(film_stack(), K0, complex(1.4302, 1e-5))
    n_slab = solve_mode(slab_stack(), K0, complex(1.4340, 0.0))
    return {
        "film_planar": n_film,
        "slab_planar": n_slab,
        "film_strip": solve_strip(K0, n_film**2, EPS_CLAD, W_AU_UM),
        "slab_strip": solve_strip(K0, n_slab**2, EPS_CLAD, W_CORE_UM),
    }


def supermodes(gap_um: float, seed=None):
    """Две планарные супермоды связанной пары; чётность - по знаку H_y."""
    stack = coupled_stack(gap_um)
    if seed is None:
        seed = (complex(1.4344, 1e-6), complex(1.4319, 3e-5))

    found: list[complex] = []
    n_clad = float(np.sqrt(EPS_CLAD).real)
    for s in seed:
        for d in (0.0, 3e-4, -3e-4, 8e-4, -8e-4, 2e-3, -2e-3):
            try:
                neff = solve_mode(stack, K0, complex(s.real + d, s.imag))
            except Exception:
                continue
            if not np.isfinite(neff.real) or neff.real <= n_clad + 1e-8:
                continue
            if neff.real > 1.4400:
                continue
            if any(abs(neff - old) < 5e-7 for old in found):
                continue
            found.append(neff)

    if len(found) < 2:
        return None
    found.sort(key=lambda z: -z.real)
    upper, lower = found[0], found[1]

    edges = stack.interfaces()
    z_probe = np.array([0.5 * (edges[0] + edges[1]), 0.5 * (edges[2] + edges[3])])
    h_upper = mode_fields(upper, stack, K0, z_probe)[0]
    same_sign = (h_upper[0].real * h_upper[1].real) > 0
    return (upper, lower) if same_sign else (lower, upper)


def gap_sweep(gaps: np.ndarray, iso: dict[str, complex]) -> list[dict]:
    rows: list[dict] = []
    seed = None
    for gap in gaps:
        pair = supermodes(float(gap), seed)
        if pair is None:
            rows.append({"gap": float(gap), "kappa": None})
            continue
        even, odd = pair
        seed = (even, odd)
        kappa = extract_kappa(LAMBDA_UM, even, odd, iso["film_planar"], iso["slab_planar"])
        device = CoupledPair(LAMBDA_UM, kappa, iso["film_strip"], iso["slab_strip"])
        best_z, best_p = device.optimal_length_um()
        rows.append({
            "gap": float(gap), "even": even, "odd": odd, "kappa": kappa,
            "lc_um": device.coupling_length_um,
            "max_transfer": device.max_transfer,
            "best_z": best_z, "best_p": best_p,
            "loss_db_cm": device.mean_loss_db_per_cm,
            "n_super_low": 0.5 * (iso["film_strip"].real + iso["slab_strip"].real)
            - device.half_split / device.k0,
        })
    return rows


def kappa_decay_rate(rows: list[dict]) -> tuple[float, float]:
    """Показатель спада связи по зазору и аналитическая оценка по хвосту поля."""
    valid = [r for r in rows if r.get("kappa")]
    s = np.array([r["gap"] for r in valid])
    k = np.array([r["kappa"] for r in valid])
    mask = k > 0
    gamma_fit = -float(np.polyfit(s[mask], np.log(k[mask]), 1)[0])
    n_mean = 0.5 * (valid[0]["even"].real + valid[0]["odd"].real)
    n_clad = float(np.sqrt(EPS_CLAD).real)
    gamma_analytic = float(K0 * np.sqrt(max(n_mean**2 - n_clad**2, 0.0)))
    return gamma_fit, gamma_analytic


def make_plots(rows: list[dict], device: CoupledPair, iso: dict[str, complex], nominal: dict) -> None:
    valid = [r for r in rows if r.get("kappa")]
    gaps = np.array([r["gap"] for r in valid])
    kappas = np.array([r["kappa"] for r in valid])
    lcs = np.array([r["lc_um"] for r in valid])

    fig, axs = plt.subplots(2, 2, figsize=(12.6, 9.0))

    ax = axs[0, 0]
    ax.semilogy(gaps, kappas, "o-", ms=3, color="#0B6E99")
    ax.set_xlabel("зазор s, мкм")
    ax.set_ylabel("коэффициент связи, мкм$^{-1}$")
    ax.set_title("Связь спадает экспоненциально с зазором")
    ax.grid(alpha=0.3, which="both")

    ax = axs[0, 1]
    ax.plot(gaps, lcs, "o-", ms=3, label="расчёт: точная связь и расстройка ЭДП")
    ax.plot(list(PAPER_LC), list(PAPER_LC.values()), "k*", ms=14, label="Park 2009, расчёт")
    ax.plot([5.0], [600.0], "rv", ms=10, label="Park 2009, эксперимент")
    ax.set_xlabel("зазор s, мкм")
    ax.set_ylabel("длина связи, мкм")
    ax.set_title("Длина связи против зазора")
    ax.set_yscale("log")
    ax.set_xlim(3.8, 10.6)
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8, frameon=False)

    ax = axs[1, 0]
    stack = coupled_stack(GAP_NOMINAL_UM)
    edges = stack.interfaces()
    z = np.linspace(-6.0, edges[-1] + 6.0, 6000)
    for neff, label, style in ((nominal["even"], "чётная", "-"), (nominal["odd"], "нечётная", "--")):
        hy = mode_fields(neff, stack, K0, z)[0]
        hy = hy / np.max(np.abs(hy))
        ax.plot(z, hy.real, style, lw=1.8, label=f"{label}, n = {neff.real:.6f}")
    ax.axvspan(edges[0], edges[1], color="goldenrod", alpha=0.7, label="Au")
    ax.axvspan(edges[2], edges[3], color="0.7", alpha=0.55, label="сердцевина")
    ax.axhline(0, color="0.4", lw=0.8)
    ax.set_xlabel("z, мкм (0 - нижняя граница золота)")
    ax.set_ylabel("H_y, норм.")
    ax.set_title(f"Планарные супермоды при зазоре {GAP_NOMINAL_UM:g} мкм")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, frameon=False)

    ax = axs[1, 1]
    z_um = np.linspace(0, 2400, 4000)
    cross, thru = device.transfer_curve(z_um)
    ax.plot(z_um, cross, lw=1.9, label="в диэлектрический волновод")
    ax.plot(z_um, thru, lw=1.4, ls="--", label="осталось в LR-SPP")
    best_z, best_p = device.optimal_length_um()
    ax.plot([best_z], [best_p], "o", color="crimson", ms=7,
            label=f"оптимум: {best_z:.0f} мкм, {best_p * 100:.0f} %")
    ax.plot([680], [PAPER_MAX_TRANSFER], "k*", ms=14, label="Park 2009: 98 % на 680 мкм")
    ax.plot([600], [0.86], "rv", ms=10, label="Park 2009: 86 % на 600 мкм")
    ax.set_xlabel("длина взаимодействия, мкм")
    ax.set_ylabel("доля мощности")
    ax.set_title("Перекачка мощности с учётом потерь")
    ax.set_ylim(0, 1.05)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, frameon=False)

    fig.suptitle("Park 2009: аналитическое воспроизведение ответвителя LR-SPP / диэлектрик", fontsize=12)
    fig.tight_layout()
    fig.savefig(OUT / "park2009_coupler.png", dpi=170)
    plt.close(fig)


def write_outputs(rows, device, iso, nominal, gamma_fit, gamma_analytic) -> list[str]:
    with (OUT / "park2009_gap_sweep.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["gap_um", "kappa_per_um", "Lc_um", "max_transfer", "best_length_um",
                    "best_transfer", "n_super_low", "n_even_planar", "n_odd_planar"])
        for r in rows:
            if not r.get("kappa"):
                continue
            w.writerow([f"{r['gap']:.3f}", f"{r['kappa']:.6e}", f"{r['lc_um']:.2f}",
                        f"{r['max_transfer']:.5f}", f"{r['best_z']:.1f}", f"{r['best_p']:.5f}",
                        f"{r['n_super_low']:.9f}", f"{r['even'].real:.9f}", f"{r['odd'].real:.9f}"])

    best_z, best_p = device.optimal_length_um()
    n_clad = float(np.sqrt(EPS_CLAD).real)
    lines: list[str] = []
    add = lines.append
    add("Park et al., Opt. Commun. 282, 4513 (2009) - аналитическое воспроизведение")
    add(f"lambda = {LAMBDA_UM} мкм; Au {W_AU_UM:g} x {T_AU_UM * 1000:.0f} нм; "
        f"сердцевина {W_CORE_UM:g} x {T_CORE_UM:g} мкм; зазор {GAP_NOMINAL_UM:g} мкм")
    add("")
    add("1. Изолированные волноводы")
    add(f"   LR-SPP, бесконечная плёнка   n_eff = {iso['film_planar'].real:.6f}"
        f"   потери {propagation_loss_db_per_cm(iso['film_planar'], LAMBDA_UM):7.2f} дБ/см")
    add(f"   LR-SPP, полоска {W_AU_UM:g} мкм       n_eff = {iso['film_strip'].real:.6f}"
        f"   потери {propagation_loss_db_per_cm(iso['film_strip'], LAMBDA_UM):7.2f} дБ/см"
        f"   (статья: {PAPER_LRSPP_LOSS:g} расчёт, 13 эксперимент)")
    add(f"   диэлектрик, бесконечный слой n_eff = {iso['slab_planar'].real:.6f}")
    add(f"   диэлектрик, полоска {W_CORE_UM:g} мкм   n_eff = {iso['slab_strip'].real:.6f}")
    add(f"   статья: оба волновода        n_eff = {PAPER_NEFF}")
    add(f"   отличие от статьи: {abs(iso['film_strip'].real - PAPER_NEFF):.5f} и "
        f"{abs(iso['slab_strip'].real - PAPER_NEFF):.5f}")
    add("")
    add("2. Связь из точного планарного расчёта")
    add(f"   при зазоре {GAP_NOMINAL_UM:g} мкм: kappa = {nominal['kappa']:.6e} мкм^-1")
    add(f"   спад по зазору, подгонка exp(-gamma s): gamma = {gamma_fit:.4f} мкм^-1")
    add(f"   аналитическая постоянная затухания в обкладке: {gamma_analytic:.4f} мкм^-1")
    add("")
    add("3. Прибор при номинальном зазоре")
    add(f"   расстройка          = {device.detuning_per_um:.6e} мкм^-1")
    add(f"   связь               = {device.kappa_per_um:.6e} мкм^-1")
    add(f"   длина связи         = {device.coupling_length_um:7.1f} мкм"
        f"   (статья: {PAPER_LC[5.0]:.0f} расчёт, 600 эксперимент)")
    add(f"   предельная перекачка= {device.max_transfer * 100:6.1f} %"
        f"   (статья: {PAPER_MAX_TRANSFER * 100:.0f} расчёт, 86 эксперимент)")
    add(f"   с учётом потерь     = {best_p * 100:6.1f} % на {best_z:.0f} мкм")
    add(f"   среднее затухание   = {device.mean_loss_db_per_cm:6.2f} дБ/см"
        f"   (статья: {PAPER_DEVICE_LOSS:g} расчёт, 7.6 эксперимент)")
    add("")
    add("4. Длина связи против зазора")
    add("   s, мкм   расчёт, мкм   статья, мкм   отличие")
    for r in rows:
        if not r.get("kappa") or r["gap"] not in PAPER_LC:
            continue
        ref = PAPER_LC[r["gap"]]
        add(f"   {r['gap']:5.1f}   {r['lc_um']:11.0f}   {ref:11.0f}   "
            f"{100 * (r['lc_um'] - ref) / ref:+6.1f} %")
    add("")
    add("5. Нижняя супермода прибора относительно обкладки")
    add(f"   показатель обкладки = {n_clad:.4f}")
    below = [r["gap"] for r in rows if r.get("kappa") and r["n_super_low"] < n_clad]
    if below:
        add(f"   уходит под обкладку при зазоре менее {max(below):.2f} мкм "
            f"(статья: {PAPER_CUTOFF_GAP:g} мкм)")
    else:
        lowest = min(r["n_super_low"] for r in rows if r.get("kappa"))
        add(f"   в диапазоне {rows[0]['gap']:.1f}-{rows[-1]['gap']:.1f} мкм остаётся выше обкладки, "
            f"минимум {lowest:.6f}")
        add(f"   статья сообщает отсечку при {PAPER_CUTOFF_GAP:g} мкм - расхождение, см. README")
    (OUT / "park2009_summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return lines


def main() -> int:
    iso = isolated_modes()
    gaps = np.round(np.arange(3.0, 12.01, 0.25), 3)
    rows = gap_sweep(gaps, iso)

    nominal = next(r for r in rows if abs(r["gap"] - GAP_NOMINAL_UM) < 1e-9)
    device = CoupledPair(LAMBDA_UM, nominal["kappa"], iso["film_strip"], iso["slab_strip"])
    gamma_fit, gamma_analytic = kappa_decay_rate(rows)

    make_plots(rows, device, iso, nominal)
    for line in write_outputs(rows, device, iso, nominal, gamma_fit, gamma_analytic):
        print(line)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
