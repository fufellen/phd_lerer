"""Аналитическое воспроизведение Park et al., Opt. Express 19, 21605 (2011).

Статья: плазмонный преобразователь размера моды. Широкий волновод ДМД (IMI,
одна плёнка золота) переходит боковым сужением в узкий волновод ДМДМД (IMIMI,
две плёнки золота с диэлектриком между ними). DOI 10.1364/OE.19.021605.

Что воспроизводится и с чем сверяется (числа статьи):
  n_eff рабочих мод ДМД и ДМДМД       ~1.450-1.452 при ширинах 2-6 мкм
  n_eff паразитной моды Sa0            1.59-1.613
  потери ДМД, полоска 6 мкм            4.38 дБ/см (измерено, cutback)
  потери ДМД, полоска 3 мкм            1.89 дБ/см (измерено)
  потери ДМДМД, полоска 3 мкм          5.01 дБ/см (измерено)
  размер моды ДМД 6 мкм (расчёт, 1/e)  11.0 x 10.5 мкм
  размер моды ДМДМД 3 мкм (1/e)        6.9 x 6.6 мкм
  связь рабочей моды с рабочей         ~0.7 дБ (расчёт)
  связь рабочей моды с паразитной      ~40 дБ (расчёт)

Важно про потери. Числа статьи по потерям - ИЗМЕРЕННЫЕ методом cutback. Park
2008 на той же технологии показал, что измеренные потери в 3-4 раза выше
расчётных на объёмных постоянных золота. Поэтому расхождение расчёта с этими
числами ожидаемо и не является ошибкой модели; сравнивать надо тенденции по
ширине и отношение потерь ДМД к ДМДМД.

Не воспроизводится здесь. Потери самого сужающегося перехода (1.16 дБ при 27
градусах и немонотонная зависимость от длины) требуют расчёта распространения
вдоль перехода - разложения по локальным модам. Это следующий этап, см. README.

Запуск:
    python lrspp_coupling/scripts/park2011_mode_size_converter.py
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
    Stack,
    materials,
    mode_fields,
    mode_width_1e,
    orthogonality_residual,
    overlap_power,
    propagation_loss_db_per_cm,
    solve_mode,
    solve_strip,
)
from slabmodes.eim import horizontal_residual_tm  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "results"
OUT.mkdir(parents=True, exist_ok=True)

LAMBDA_UM = 1.55
K0 = materials.k0_from_lambda(LAMBDA_UM)

EPS_CLAD = materials.EPS_ZPU450  # ZPU450, n = 1.450
EPS_AU = materials.EPS_AU_1550
T_AU_UM = 0.014  # толщина каждой плёнки золота
T_GAP_UM = 0.500  # центральный диэлектрик ДМДМД

W_IMI_UM = 6.0
W_IMIMI_UM = 3.0
WIDTHS = np.array([2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0])

PAPER_PL_IMI = {3.0: 1.89, 6.0: 4.38}  # дБ/см, измерено
PAPER_PL_IMIMI = {3.0: 5.01}
PAPER_CL_PMF = {3.0: 4.52, 6.0: 0.77}  # дБ на грань, измерено
PAPER_SIZE_IMI = (11.0, 10.5)  # мкм, расчёт, уровень 1/e
PAPER_SIZE_IMIMI = (6.9, 6.6)
PAPER_OVERLAP_WORK_DB = 0.7
PAPER_OVERLAP_PARASITIC_DB = 40.0


def imi_stack() -> Stack:
    """ДМД: полимер | Au | полимер."""
    return Stack(eps=(EPS_CLAD, EPS_AU, EPS_CLAD), thickness=(T_AU_UM,),
                 names=("полимер", "Au", "полимер"))


def imimi_stack() -> Stack:
    """ДМДМД: полимер | Au | полимер | Au | полимер."""
    return Stack(
        eps=(EPS_CLAD, EPS_AU, EPS_CLAD, EPS_AU, EPS_CLAD),
        thickness=(T_AU_UM, T_GAP_UM, T_AU_UM),
        names=("полимер", "Au", "зазор", "Au", "полимер"),
    )


def parity_of(neff: complex, stack: Stack) -> str:
    """Чётность H_y относительно центра структуры."""
    edges = stack.interfaces()
    centre = 0.5 * (edges[0] + edges[-1])
    half = 0.35 * (edges[-1] - edges[0]) + 0.2
    z = np.array([centre - half, centre + half])
    h = mode_fields(neff, stack, K0, z)[0]
    return "чётная" if (h[0].real * h[1].real) > 0 else "нечётная"


def imimi_modes() -> dict[str, complex]:
    """Рабочая Ss0 и паразитная Sa0 моды ДМДМД."""
    stack = imimi_stack()
    n_clad = float(np.sqrt(EPS_CLAD).real)
    found: list[complex] = []
    guesses = [complex(1.4505, 1e-5), complex(1.4515, 3e-5), complex(1.4530, 1e-4),
               complex(1.470, 1e-3), complex(1.520, 3e-3), complex(1.560, 5e-3),
               complex(1.590, 8e-3), complex(1.610, 1e-2), complex(1.650, 2e-2),
               complex(1.700, 3e-2)]
    for guess in guesses:
        try:
            neff = solve_mode(stack, K0, guess)
        except Exception:
            continue
        if not np.isfinite(neff.real) or neff.real <= n_clad + 1e-7 or neff.real > 3.0:
            continue
        if any(abs(neff - old) < 1e-6 for old in found):
            continue
        found.append(neff)
    found.sort(key=lambda z: z.real)
    out: dict[str, complex] = {}
    for neff in found:
        out[f"{parity_of(neff, stack)} n={neff.real:.6f}"] = neff
    return out


def horizontal_profile(neff_strip: complex, eps_core: complex, width_um: float,
                       x: np.ndarray) -> np.ndarray:
    """Поперечный профиль поля по горизонтальной задаче ЭДП (чётное решение)."""
    half = width_um / 2.0
    u = K0 * np.sqrt(eps_core - neff_strip * neff_strip + 0j)
    alpha = K0 * np.sqrt(neff_strip * neff_strip - EPS_CLAD + 0j)
    if alpha.real < 0:
        alpha = -alpha
    inside = np.abs(x) <= half
    prof = np.zeros_like(x, dtype=complex)
    prof[inside] = np.cos(u * x[inside])
    edge = np.cos(u * half)
    prof[~inside] = edge * np.exp(-alpha * (np.abs(x[~inside]) - half))
    return prof


def mode_size_2d(neff_planar: complex, stack: Stack, width_um: float) -> tuple[float, float]:
    """Размеры моды по уровню 1/e: вертикальный из точного профиля, горизонтальный из ЭДП."""
    edges = stack.interfaces()
    z = np.linspace(edges[0] - 40.0, edges[-1] + 40.0, 200001)
    hy = mode_fields(neff_planar, stack, K0, z)[0]
    vertical = mode_width_1e(hy, z)

    n_strip = solve_strip(K0, neff_planar**2, EPS_CLAD, width_um)
    x = np.linspace(-40.0, 40.0, 200001)
    prof = horizontal_profile(n_strip, neff_planar**2, width_um, x)
    horizontal = mode_width_1e(prof, x)
    return horizontal, vertical


def width_sweep() -> list[dict]:
    n_imi = solve_mode(imi_stack(), K0, complex(1.4512, 2e-5))
    modes = imimi_modes()
    work = min((v for v in modes.values()), key=lambda z: z.real)

    rows: list[dict] = []
    for w in WIDTHS:
        row: dict = {"width": float(w)}
        try:
            n_s = solve_strip(K0, n_imi**2, EPS_CLAD, float(w))
            row["imi"] = n_s
            row["imi_loss"] = propagation_loss_db_per_cm(n_s, LAMBDA_UM)
        except RuntimeError:
            row["imi"] = None
        try:
            n_s = solve_strip(K0, work**2, EPS_CLAD, float(w))
            row["imimi"] = n_s
            row["imimi_loss"] = propagation_loss_db_per_cm(n_s, LAMBDA_UM)
        except RuntimeError:
            row["imimi"] = None
        rows.append(row)
    return rows


def vertical_overlaps() -> dict[str, float]:
    """Перекрытие рабочей моды ДМД с модами ДМДМД по вертикальному сечению."""
    st_imi, st_imimi = imi_stack(), imimi_stack()
    n_imi = solve_mode(st_imi, K0, complex(1.4512, 2e-5))
    modes = imimi_modes()

    # общая сетка с центрами обеих структур в нуле
    z = np.linspace(-60.0, 60.0, 400001)
    c_imi = 0.5 * (st_imi.interfaces()[0] + st_imi.interfaces()[-1])
    c_imm = 0.5 * (st_imimi.interfaces()[0] + st_imimi.interfaces()[-1])

    h_imi = mode_fields(n_imi, st_imi, K0, z + c_imi)[0]
    eps_imi = st_imi.eps_at(z + c_imi)

    out: dict[str, float] = {}
    for label, neff in modes.items():
        h = mode_fields(neff, st_imimi, K0, z + c_imm)[0]
        eps = st_imimi.eps_at(z + c_imm)
        eta = overlap_power((h_imi, n_imi, eps_imi), (h, neff, eps), z)
        out[label] = eta
    return out


def make_plots(rows: list[dict], modes: dict[str, complex], overlaps: dict[str, float]) -> None:
    st_imi, st_imimi = imi_stack(), imimi_stack()
    n_imi = solve_mode(st_imi, K0, complex(1.4512, 2e-5))

    fig, axs = plt.subplots(2, 2, figsize=(12.6, 9.0))

    ax = axs[0, 0]
    z = np.linspace(-8.0, 8.0, 8000)
    h = mode_fields(n_imi, st_imi, K0, z + 0.5 * st_imi.interfaces()[-1])[0]
    ax.plot(z, (h / np.max(np.abs(h))).real, lw=2.0, label=f"ДМД, n = {n_imi.real:.6f}")
    c = 0.5 * st_imimi.interfaces()[-1]
    for label, neff in sorted(modes.items(), key=lambda kv: kv[1].real):
        hh = mode_fields(neff, st_imimi, K0, z + c)[0]
        ax.plot(z, (hh / np.max(np.abs(hh))).real, lw=1.5, ls="--", label=f"ДМДМД, {label}")
    ax.axhline(0, color="0.4", lw=0.8)
    ax.set_xlim(-4, 4)
    ax.set_xlabel("z, мкм (0 - центр структуры)")
    ax.set_ylabel("H_y, норм.")
    ax.set_title("Вертикальные профили рабочих и паразитных мод")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, frameon=False)

    ax = axs[0, 1]
    widths = [r["width"] for r in rows if r.get("imi")]
    imi_n = [r["imi"].real for r in rows if r.get("imi")]
    ax.plot(widths, imi_n, "o-", ms=4, label="ДМД, полоска")
    w2 = [r["width"] for r in rows if r.get("imimi")]
    if w2:
        ax.plot(w2, [r["imimi"].real for r in rows if r.get("imimi")], "s-", ms=4, label="ДМДМД, полоска")
    ax.axhline(np.sqrt(EPS_CLAD).real, color="0.3", ls=":", label="показатель обкладки")
    ax.set_xlabel("ширина полоски, мкм")
    ax.set_ylabel("Re n_eff")
    ax.set_title("Эффективный показатель против ширины")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, frameon=False)

    ax = axs[1, 0]
    ax.plot(widths, [r["imi_loss"] for r in rows if r.get("imi")], "o-", ms=4, label="ДМД, расчёт")
    if w2:
        ax.plot(w2, [r["imimi_loss"] for r in rows if r.get("imimi")], "s-", ms=4, label="ДМДМД, расчёт")
    ax.plot(list(PAPER_PL_IMI), list(PAPER_PL_IMI.values()), "k*", ms=14, label="ДМД, измерено")
    ax.plot(list(PAPER_PL_IMIMI), list(PAPER_PL_IMIMI.values()), "rv", ms=10, label="ДМДМД, измерено")
    ax.set_xlabel("ширина полоски, мкм")
    ax.set_ylabel("потери, дБ/см")
    ax.set_title("Потери распространения: расчёт против измерений")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, frameon=False)

    ax = axs[1, 1]
    labels = list(overlaps)
    values = [max(overlaps[k], 1e-16) for k in labels]
    db = [-10.0 * np.log10(v) for v in values]
    colors = ["#0B6E99" if d < 10 else "#B23A48" for d in db]
    ax.barh(range(len(labels)), db, color=colors)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels([f"ДМД -> {k}" for k in labels], fontsize=7)
    ax.axvline(PAPER_OVERLAP_WORK_DB, color="k", ls="--", lw=1.2, label="статья: 0.7 дБ")
    ax.axvline(PAPER_OVERLAP_PARASITIC_DB, color="crimson", ls=":", lw=1.2, label="статья: 40 дБ")
    ax.set_xlabel("потери на стыке, дБ")
    ax.set_title("Связь моды ДМД с модами ДМДМД")
    ax.grid(alpha=0.3, axis="x")
    ax.legend(fontsize=8, frameon=False)

    fig.suptitle("Park 2011: аналитическое воспроизведение преобразователя размера моды", fontsize=12)
    fig.tight_layout()
    fig.savefig(OUT / "park2011_converter.png", dpi=170)
    plt.close(fig)


def main() -> int:
    st_imi, st_imimi = imi_stack(), imimi_stack()
    n_imi = solve_mode(st_imi, K0, complex(1.4512, 2e-5))
    modes = imimi_modes()
    rows = width_sweep()
    overlaps = vertical_overlaps()

    size_imi = mode_size_2d(n_imi, st_imi, W_IMI_UM)
    work = min(modes.values(), key=lambda z: z.real)
    size_imimi = mode_size_2d(work, st_imimi, W_IMIMI_UM)

    lines: list[str] = []
    add = lines.append
    add("Park et al., Opt. Express 19, 21605 (2011) - аналитическое воспроизведение")
    add(f"lambda = {LAMBDA_UM} мкм; Au {T_AU_UM * 1000:.0f} нм; "
        f"центральный слой ДМДМД {T_GAP_UM * 1000:.0f} нм; полимер n = {np.sqrt(EPS_CLAD).real:.3f}")
    add("")
    add("1. Планарные моды")
    add(f"   ДМД, рабочая мода     n_eff = {n_imi.real:.6f}{n_imi.imag:+.3e}i   "
        f"потери {propagation_loss_db_per_cm(n_imi, LAMBDA_UM):7.3f} дБ/см   чётность {parity_of(n_imi, st_imi)}")
    for label, neff in sorted(modes.items(), key=lambda kv: kv[1].real):
        add(f"   ДМДМД {label:26s} потери {propagation_loss_db_per_cm(neff, LAMBDA_UM):8.3f} дБ/см")
    add(f"   статья: рабочие моды ~1.450-1.452, паразитная Sa0 1.59-1.613")
    add("")
    add("2. Перекрытие рабочей моды ДМД с модами ДМДМД (вертикальное сечение)")
    for label, eta in sorted(overlaps.items(), key=lambda kv: -kv[1]):
        db = -10.0 * np.log10(max(eta, 1e-16))
        add(f"   -> {label:30s} eta = {eta:.6e}   {db:8.2f} дБ")
    add(f"   статья: рабочая ~{PAPER_OVERLAP_WORK_DB} дБ, паразитная ~{PAPER_OVERLAP_PARASITIC_DB} дБ")
    add("")
    add("3. Размеры моды по уровню 1/e (горизонталь x вертикаль), мкм")
    add(f"   ДМД, полоска {W_IMI_UM:g} мкм    расчёт {size_imi[0]:5.2f} x {size_imi[1]:5.2f}"
        f"   статья {PAPER_SIZE_IMI[0]} x {PAPER_SIZE_IMI[1]}")
    add(f"   ДМДМД, полоска {W_IMIMI_UM:g} мкм  расчёт {size_imimi[0]:5.2f} x {size_imimi[1]:5.2f}"
        f"   статья {PAPER_SIZE_IMIMI[0]} x {PAPER_SIZE_IMIMI[1]}")
    area_ratio = (size_imimi[0] * size_imimi[1]) / (size_imi[0] * size_imi[1])
    paper_ratio = (PAPER_SIZE_IMIMI[0] * PAPER_SIZE_IMIMI[1]) / (PAPER_SIZE_IMI[0] * PAPER_SIZE_IMI[1])
    add(f"   сжатие по площади: расчёт {area_ratio * 100:.1f} %, статья {paper_ratio * 100:.1f} %")
    add("")
    add("4. Потери против ширины полоски, дБ/см")
    add("   W, мкм   ДМД расчёт   ДМД измерено   ДМДМД расчёт   ДМДМД измерено")
    for r in rows:
        imi_l = f"{r['imi_loss']:10.3f}" if r.get("imi") else "         -"
        imm_l = f"{r['imimi_loss']:12.3f}" if r.get("imimi") else "           -"
        ref_i = f"{PAPER_PL_IMI[r['width']]:12.2f}" if r["width"] in PAPER_PL_IMI else "           -"
        ref_m = f"{PAPER_PL_IMIMI[r['width']]:14.2f}" if r["width"] in PAPER_PL_IMIMI else "             -"
        add(f"   {r['width']:5.1f}   {imi_l}   {ref_i}   {imm_l}   {ref_m}")
    add("")
    add("   Измеренные значения получены методом cutback. Park 2008 на той же")
    add("   технологии показал превышение измеренных потерь над расчётными в 3-4 раза,")
    add("   поэтому прямое совпадение здесь не ожидается: сравниваются тенденции.")

    with (OUT / "park2011_width_sweep.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["width_um", "imi_neff_real", "imi_neff_imag", "imi_loss_db_cm",
                    "imimi_neff_real", "imimi_neff_imag", "imimi_loss_db_cm"])
        for r in rows:
            w.writerow([
                f"{r['width']:.2f}",
                f"{r['imi'].real:.9f}" if r.get("imi") else "",
                f"{r['imi'].imag:.6e}" if r.get("imi") else "",
                f"{r['imi_loss']:.4f}" if r.get("imi") else "",
                f"{r['imimi'].real:.9f}" if r.get("imimi") else "",
                f"{r['imimi'].imag:.6e}" if r.get("imimi") else "",
                f"{r['imimi_loss']:.4f}" if r.get("imimi") else "",
            ])
    (OUT / "park2011_summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

    make_plots(rows, modes, overlaps)
    for line in lines:
        print(line)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
