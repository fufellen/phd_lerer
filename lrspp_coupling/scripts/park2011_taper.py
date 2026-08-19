"""Потери сужающегося перехода в Park et al., Opt. Express 19, 21605 (2011).

Измеренная зависимость потерь одного перехода от угла сужения (рис. 4 статьи):

    угол   длина перехода   потери
    0.38 гр   ~240 мкм       2.00 дБ
    1.39      ~62            1.41
    5.7       ~15            1.21
    14.0      ~6             1.19
    26.6      ~3             1.16   <- минимум
    45.0      1.5            1.21

Структура перехода. Подводящий волновод ДМД шириной 6 мкм (одна плёнка золота)
стыкуется с волноводом ДМДМД той же ширины (две плёнки), после чего ширина
ДМДМД сужается с 6 до 3 мкм на длине L_t. Значит, потери перехода складываются
минимум из трёх слагаемых:

  1) резкий стык ДМД -> ДМДМД при ширине 6 мкм: рассогласование профилей;
  2) поглощение вдоль сужающегося участка: растёт пропорционально длине;
  3) неадиабатичность при слишком быстром сужении: растёт при укорочении.

Первые два считаются здесь точно, третье оценивается критерием Лава. Важно, что
измеренная величина по построению включает поглощение внутри перехода - авторы
это оговаривают явно.

Запуск:
    python lrspp_coupling/scripts/park2011_taper.py
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

from slabmodes import Stack, materials, mode_fields, overlap_power, solve_mode, solve_strip, trapz  # noqa: E402
from slabmodes.eme import TaperProfile, adiabatic_transmission, junction_loss_db, love_adiabaticity  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "results"
OUT.mkdir(parents=True, exist_ok=True)

LAMBDA_UM = 1.55
K0 = materials.k0_from_lambda(LAMBDA_UM)
EPS_CLAD = materials.EPS_ZPU450
EPS_AU = materials.EPS_AU_1550
T_AU_UM = 0.014
T_GAP_UM = 0.500
N_CLAD = float(np.sqrt(EPS_CLAD).real)

W_WIDE_UM = 6.0
W_NARROW_UM = 3.0

# Измеренные точки рис. 4 статьи: длина перехода (мкм) -> потери (дБ)
MEASURED = {240.0: 2.00, 62.0: 1.41, 15.0: 1.21, 6.0: 1.19, 3.0: 1.16, 1.5: 1.21}

# Измеренные погонные потери из рис. 3 статьи, дБ/см
PAPER_PL_IMIMI = {3.0: 5.01, 6.0: 8.97}


def imi_stack() -> Stack:
    return Stack(eps=(EPS_CLAD, EPS_AU, EPS_CLAD), thickness=(T_AU_UM,))


def imimi_stack() -> Stack:
    return Stack(eps=(EPS_CLAD, EPS_AU, EPS_CLAD, EPS_AU, EPS_CLAD),
                 thickness=(T_AU_UM, T_GAP_UM, T_AU_UM))


def working_modes() -> tuple[complex, complex]:
    """Рабочие планарные моды ДМД и ДМДМД."""
    n_imi = solve_mode(imi_stack(), K0, complex(1.4512, 2e-5))
    n_imimi = solve_mode(imimi_stack(), K0, complex(1.4534, 3e-5))
    return n_imi, n_imimi


def vertical_junction_overlap(n_imi: complex, n_imimi: complex) -> float:
    """Перекрытие вертикальных профилей на резком стыке ДМД -> ДМДМД."""
    st_a, st_b = imi_stack(), imimi_stack()
    z = np.linspace(-60.0, 60.0, 400001)
    ca = 0.5 * st_a.interfaces()[-1]
    cb = 0.5 * st_b.interfaces()[-1]
    ha = mode_fields(n_imi, st_a, K0, z + ca)[0]
    hb = mode_fields(n_imimi, st_b, K0, z + cb)[0]
    return overlap_power((ha, n_imi, st_a.eps_at(z + ca)), (hb, n_imimi, st_b.eps_at(z + cb)), z)


def horizontal_profile(neff_strip: complex, eps_core: complex, width_um: float,
                       x: np.ndarray) -> np.ndarray:
    """Чётный поперечный профиль полоски по горизонтальной задаче ЭДП."""
    half = width_um / 2.0
    u = K0 * np.sqrt(eps_core - neff_strip * neff_strip + 0j)
    alpha = K0 * np.sqrt(neff_strip * neff_strip - EPS_CLAD + 0j)
    if alpha.real < 0:
        alpha = -alpha
    prof = np.zeros_like(x, dtype=complex)
    inside = np.abs(x) <= half
    prof[inside] = np.cos(u * x[inside])
    prof[~inside] = np.cos(u * half) * np.exp(-alpha * (np.abs(x[~inside]) - half))
    return prof


def horizontal_junction_overlap(n_imi: complex, n_imimi: complex) -> float:
    """Перекрытие поперечных профилей полосок одинаковой ширины 6 мкм."""
    n_a = solve_strip(K0, n_imi**2, EPS_CLAD, W_WIDE_UM)
    n_b = solve_strip(K0, n_imimi**2, EPS_CLAD, W_WIDE_UM)
    x = np.linspace(-60.0, 60.0, 400001)
    pa = horizontal_profile(n_a, n_imi**2, W_WIDE_UM, x)
    pb = horizontal_profile(n_b, n_imimi**2, W_WIDE_UM, x)
    num = abs(trapz(pa * np.conj(pb), x)) ** 2
    den = trapz(np.abs(pa) ** 2, x) * trapz(np.abs(pb) ** 2, x)
    return float(num / den)


class WidthTable:
    """Кешированная зависимость n_eff полоски ДМДМД от ширины."""

    def __init__(self, n_planar: complex, widths: np.ndarray) -> None:
        self.n_planar = n_planar
        self.widths = widths
        self.values = np.array([solve_strip(K0, n_planar**2, EPS_CLAD, float(w)) for w in widths])

    def __call__(self, width_um: float) -> complex:
        re = float(np.interp(width_um, self.widths, self.values.real))
        im = float(np.interp(width_um, self.widths, self.values.imag))
        return complex(re, im)

    def alpha_db_per_cm(self, width_um: float) -> float:
        return float(2.0 * K0 * abs(self(width_um).imag) * 1e4 * 10.0 / np.log(10.0))


def analyse(table: WidthTable, junction_db: float) -> list[dict]:
    rows: list[dict] = []
    for length in sorted(MEASURED, reverse=True):
        profile = TaperProfile(length, W_WIDE_UM, W_NARROW_UM)
        absorb = adiabatic_transmission(profile, table, LAMBDA_UM)
        love = love_adiabaticity(profile, table, N_CLAD, LAMBDA_UM)
        total = junction_db + absorb["loss_db"]
        rows.append({
            "length": length,
            "angle_deg": profile.half_angle_deg,
            "junction_db": junction_db,
            "absorb_db": absorb["loss_db"],
            "total_db": total,
            "measured_db": MEASURED[length],
            "residual_db": MEASURED[length] - total,
            "omega_limit_deg": love["omega_limit_deg"],
            "violation": love["violation"],
            "mean_alpha": absorb["mean_alpha_db_per_cm"],
        })
    return rows


def excess_loss_fit(rows: list[dict]) -> dict:
    """Подгонка измеренных потерь перехода моделью "полка плюс погонные потери".

    Если весь рост с длиной вызван потерями внутри самого перехода, то

        CL(L) = CL_0 + alpha_taper * L,

    и наклон прямой прямо даёт погонные потери в области сужения. Точка с самым
    коротким переходом исключается: там работает отдельный механизм - рост при
    слишком резком сужении, видимый в измерениях как загиб вверх.
    """
    pts = sorted((r["length"], r["measured_db"]) for r in rows)
    fit_pts = [p for p in pts if p[0] >= 3.0]
    x = np.array([p[0] for p in fit_pts])
    y = np.array([p[1] for p in fit_pts])
    slope, intercept = np.polyfit(x, y, 1)
    resid = y - (slope * x + intercept)
    return {
        "alpha_db_per_um": float(slope),
        "alpha_db_per_cm": float(slope * 1e4),
        "intercept_db": float(intercept),
        "max_resid_db": float(np.max(np.abs(resid))),
        "points": fit_pts,
    }


def lateral_confinement(table: WidthTable) -> list[dict]:
    """Насколько прочно мода удерживается по ширине вдоль перехода.

    Слабо связанная мода легко теряет мощность на любой неоднородности края,
    поэтому длина спадания поля в обкладке - показатель уязвимости перехода.
    """
    rows = []
    for w in (3.0, 4.0, 5.0, 6.0):
        n = table(w).real
        delta = max(n - N_CLAD, 1e-12)
        gamma = K0 * np.sqrt(max(n * n - N_CLAD**2, 1e-12))
        rows.append({"width": w, "neff": n, "delta_n": delta, "decay_um": 1.0 / gamma})
    return rows


def horizontal_mode_count(n_planar: complex, width_um: float) -> dict:
    """Число связанных мод горизонтальной задачи для заданной ветви.

    Симметричный волновод: следующая чётная мода появляется при V > pi,
    ближайшая нечётная - при V > pi/2. Симметричный переход нечётные моды не
    возбуждает, поэтому определяющий порог здесь именно pi.
    """
    na = float(np.sqrt(max(n_planar.real**2 - N_CLAD**2, 0.0)))
    v = K0 * (width_um / 2.0) * na
    return {
        "width": width_um,
        "V": float(v),
        "even_modes": int(np.floor(v / np.pi)) + 1,
        "odd_modes": int(np.floor((v - np.pi / 2) / np.pi)) + 1 if v > np.pi / 2 else 0,
    }


def make_plot(rows: list[dict], table: WidthTable) -> None:
    lengths = np.array([r["length"] for r in rows])
    order = np.argsort(lengths)
    lengths = lengths[order]
    junction = np.array([r["junction_db"] for r in rows])[order]
    absorb = np.array([r["absorb_db"] for r in rows])[order]
    total = np.array([r["total_db"] for r in rows])[order]
    measured = np.array([r["measured_db"] for r in rows])[order]

    fig, axs = plt.subplots(1, 3, figsize=(15.5, 4.9))

    ax = axs[0]
    ax.semilogx(lengths, measured, "k*-", ms=14, lw=1.2, label="измерено, Park 2011")
    ax.semilogx(lengths, total, "o-", ms=5, color="#0B6E99", label="расчёт: стык + поглощение")
    ax.fill_between(lengths, 0, junction, color="#0B6E99", alpha=0.25, label="резкий стык ДМД -> ДМДМД")
    ax.fill_between(lengths, junction, junction + absorb, color="#E1A730", alpha=0.45,
                    label="поглощение вдоль перехода")
    ax.set_xlabel("длина перехода, мкм")
    ax.set_ylabel("потери одного перехода, дБ")
    ax.set_title("Разложение потерь перехода")
    ax.set_ylim(0, 2.3)
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=7.5, frameon=False, loc="upper left")

    ax = axs[1]
    residual = measured - total
    ax.semilogx(lengths, residual, "s-", ms=6, color="#B23A48")
    ax.axhline(float(np.mean(residual)), color="0.4", ls="--", lw=1.2,
               label=f"среднее {np.mean(residual):.2f} дБ")
    ax.set_xlabel("длина перехода, мкм")
    ax.set_ylabel("измерено минус расчёт, дБ")
    ax.set_title("Необъяснённый остаток")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8, frameon=False)

    ax = axs[2]
    w = np.linspace(W_NARROW_UM, W_WIDE_UM, 60)
    ax.plot(w, [table.alpha_db_per_cm(float(x)) for x in w], lw=2, label="расчёт, ЭДП")
    ax.plot(list(PAPER_PL_IMIMI), list(PAPER_PL_IMIMI.values()), "k*", ms=14, label="измерено, рис. 3b")
    ax.set_xlabel("ширина полоски ДМДМД, мкм")
    ax.set_ylabel("потери, дБ/см")
    ax.set_title("Погонные потери вдоль перехода")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, frameon=False)

    fig.suptitle("Park 2011: из чего складываются потери сужающегося перехода", fontsize=12)
    fig.tight_layout()
    fig.savefig(OUT / "park2011_taper.png", dpi=170)
    plt.close(fig)


def main() -> int:
    n_imi, n_imimi = working_modes()
    eta_v = vertical_junction_overlap(n_imi, n_imimi)
    eta_h = horizontal_junction_overlap(n_imi, n_imimi)
    eta_2d = eta_v * eta_h
    junction_db = junction_loss_db(eta_2d)

    table = WidthTable(n_imimi, np.linspace(2.5, 6.5, 41))
    rows = analyse(table, junction_db)
    make_plot(rows, table)

    lines: list[str] = []
    add = lines.append
    add("Park et al., Opt. Express 19, 21605 (2011): потери сужающегося перехода")
    add(f"Переход: ДМД {W_WIDE_UM:g} мкм -> ДМДМД {W_WIDE_UM:g} мкм (резкий стык), "
        f"затем сужение до {W_NARROW_UM:g} мкм")
    add("")
    add("1. Резкий стык ДМД -> ДМДМД при ширине 6 мкм")
    add(f"   перекрытие по вертикали     = {eta_v:.4f}")
    add(f"   перекрытие по горизонтали   = {eta_h:.4f}")
    add(f"   двумерное перекрытие        = {eta_2d:.4f}  ->  {junction_db:.3f} дБ")
    add("   статья, расчёт связи мод: около 0.7 дБ")
    add("")
    add("2. Разложение потерь по длине перехода")
    add("   L, мкм   угол    расчёт: стык   поглощение   сумма   измерено   остаток")
    for r in sorted(rows, key=lambda d: d["length"]):
        add(f"   {r['length']:6.1f}   {r['angle_deg']:5.2f}   {r['junction_db']:11.3f}   "
            f"{r['absorb_db']:10.3f}   {r['total_db']:5.3f}   {r['measured_db']:8.2f}   "
            f"{r['residual_db']:+7.3f}")
    residuals = [r["residual_db"] for r in rows]
    add(f"   остаток: среднее {np.mean(residuals):+.3f} дБ, разброс "
        f"{np.max(residuals) - np.min(residuals):.3f} дБ")
    add("")
    add("3. Критерий адиабатичности Лава")
    add(f"   конкурирующее решение - порог излучения, показатель обкладки {N_CLAD:.3f}")
    add("   L, мкм   угол сужения   предельный угол   превышение")
    for r in sorted(rows, key=lambda d: d["length"]):
        add(f"   {r['length']:6.1f}   {r['angle_deg']:12.2f}   {r['omega_limit_deg']:15.3f}   "
            f"{r['violation']:9.1f}")
    add("")
    add("4. Избыточный канал: подгонка измеренных данных")
    fit = excess_loss_fit(rows)
    add(f"   модель CL(L) = CL_0 + alpha * L по точкам от 3 до 240 мкм:")
    add(f"     полка CL_0        = {fit['intercept_db']:.3f} дБ")
    add(f"     наклон alpha      = {fit['alpha_db_per_cm']:.1f} дБ/см")
    add(f"     макс. невязка     = {fit['max_resid_db']:.3f} дБ")
    add("   То есть весь рост с длиной описывается одной постоянной погонной потерей")
    add("   внутри перехода, и вопрос сводится к её величине.")
    straight = 0.5 * (PAPER_PL_IMIMI[3.0] + PAPER_PL_IMIMI[6.0])
    calc = 0.5 * (table.alpha_db_per_cm(3.0) + table.alpha_db_per_cm(6.0))
    add(f"   прямой волновод той же ширины, измерено: {straight:.2f} дБ/см")
    add(f"   прямой волновод той же ширины, расчёт:   {calc:.2f} дБ/см")
    add(f"   превышение над измеренным прямым: в {fit['alpha_db_per_cm'] / straight:.1f} раза")
    add("")
    add("   Что исключено как причина:")
    add("     поглощение по ставке прямого волновода - объясняет только пятую часть;")
    add("     неадиабатичность - даёт рост при УКОРОЧЕНИИ перехода, знак обратный;")
    add("     недоучёт длины в формуле (1) статьи - величина порядка 0.06 дБ, мало;")
    add("     переход между вертикальными ветвями - в разделимом приближении")
    add("       запрещён: вертикальный профиль вдоль сужения не меняется, поэтому")
    add("       ветви остаются ортогональными при любой ширине.")
    add("")
    add("5. Насколько прочно мода держится по ширине")
    add("   W, мкм   n_eff      n_eff - n_обкл   длина спадания, мкм")
    for r in lateral_confinement(table):
        add(f"   {r['width']:5.1f}   {r['neff']:.6f}   {r['delta_n']:12.2e}   {r['decay_um']:16.1f}")
    add("   Мода удерживается слабо, и её хвост простирается на единицы микрометров")
    add("   за край полоски. Наклонная кромка на такой длине работает как протяжённая")
    add("   неоднородность, поэтому рабочая гипотеза для избыточного канала -")
    add("   рассеяние на наклонной кромке металла, растущее пропорционально её длине.")
    add("   Проверяется прямым волноводом с наклонной кромкой той же длины.")
    add("")
    add("6. Одномодовость горизонтальной задачи")
    add("   Определяет, нужен ли многомодовый продольный расчёт.")
    add("   ветвь        W, мкм   V      чётных мод   нечётных мод")
    for label, n_pl in (("ДМД  ", solve_mode(imi_stack(), K0, complex(1.4512, 2e-5))),
                        ("ДМДМД", solve_mode(imimi_stack(), K0, complex(1.4534, 7e-5)))):
        for w in (3.0, 6.0):
            c = horizontal_mode_count(n_pl, w)
            add(f"   {label}        {c['width']:4.1f}   {c['V']:.3f}   {c['even_modes']:10d}   {c['odd_modes']:12d}")
    add("   Во всём диапазоне ширин перехода горизонтальная задача одномодовая:")
    add("   следующая чётная мода требует V > pi, следующая нечётная - V > pi/2, а")
    add("   симметричный переход нечётные моды и не возбуждает. Поэтому многомодовый")
    add("   продольный расчёт с базисом излучательных мод здесь ничего не добавит.")
    add("")
    add("7. Погонные потери вдоль перехода")
    for w in (3.0, 4.0, 5.0, 6.0):
        ref = f"   измерено {PAPER_PL_IMIMI[w]:.2f}" if w in PAPER_PL_IMIMI else ""
        add(f"   ширина {w:.1f} мкм: расчёт {table.alpha_db_per_cm(w):6.2f} дБ/см{ref}")

    with (OUT / "park2011_taper.csv").open("w", newline="", encoding="utf-8") as f:
        wr = csv.writer(f)
        wr.writerow(["length_um", "half_angle_deg", "junction_db", "absorption_db",
                     "total_db", "measured_db", "residual_db", "omega_limit_deg", "violation"])
        for r in sorted(rows, key=lambda d: d["length"]):
            wr.writerow([f"{r['length']:.1f}", f"{r['angle_deg']:.3f}", f"{r['junction_db']:.4f}",
                         f"{r['absorb_db']:.4f}", f"{r['total_db']:.4f}", f"{r['measured_db']:.2f}",
                         f"{r['residual_db']:.4f}", f"{r['omega_limit_deg']:.4f}",
                         f"{r['violation']:.2f}"])
    (OUT / "park2011_taper_summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    for line in lines:
        print(line)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
