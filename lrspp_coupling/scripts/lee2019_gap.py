"""Аналитическое воспроизведение Lee & Lee, J. Nanosci. Nanotechnol. 19, 6192 (2019).

Статья: LR-SPP-волновод с разрывом длиной 8 мкм и короткими обратными тейперами
по обе стороны разрыва. DOI 10.1166/jnn.2019.17010.

Что воспроизводится (измерено, lambda = 1.55 мкм):
  избыточные потери участка с разрывом 8 мкм без тейперов   1.03 дБ
  то же с тейперами до 6 мкм на длине 3 мкм                 0.76 дБ
  минимум по ширине тейпера (5 мкм)                         0.74 дБ
  выигрыш от пары тейперов                                  0.27 дБ

Физика. В разрыве металла нет, поэтому мода перестаёт быть связанной и свободно
дифрагирует в однородном полимере, а на другой стороне часть поля захватывается
обратно. Значит, потери разрыва - это перекрытие продифрагировавшего поля с
модой выходного волновода. Расчёт ведётся методом углового спектра.

Обратный тейпер расширяет моду перед разрывом. Более широкий пучок дифрагирует
медленнее (дифракционная расходимость обратно пропорциональна ширине), поэтому
теряется меньше - именно так авторы и объясняют выигрыш.

Запуск:
    python lrspp_coupling/scripts/lee2019_gap.py
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
from slabmodes.coupling import coupling_loss_db, overlap_2d, propagate_angular_spectrum, strip_mode_field  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "results"
OUT.mkdir(parents=True, exist_ok=True)

LAMBDA_UM = 1.55
K0 = materials.k0_from_lambda(LAMBDA_UM)
EPS_CLAD = materials.EPS_ZPU450
EPS_AU = materials.EPS_AU_1550
N_CLAD = float(np.sqrt(EPS_CLAD).real)

T_AU_UM = 0.020
W_GUIDE_UM = 2.0
GAP_UM = 8.0
TAPER_WIDTHS = np.array([2.0, 3.0, 4.0, 5.0, 6.0, 7.0])

MEASURED = {2.0: 1.034, 3.0: 0.804, 4.0: 0.765, 5.0: 0.744, 6.0: 0.765, 7.0: 0.814}
PAPER_PROP_DB_CM = 8.73


def stack() -> Stack:
    return Stack(eps=(EPS_CLAD, EPS_AU, EPS_CLAD), thickness=(T_AU_UM,))


def mode_field(width_um: float, x: np.ndarray, z: np.ndarray, n_planar: complex) -> np.ndarray:
    n_strip = solve_strip(K0, n_planar**2, EPS_CLAD, width_um)
    return strip_mode_field(n_planar, n_strip, stack(), K0, width_um, EPS_CLAD, x, z), n_strip


def gap_loss(width_at_gap_um: float, n_planar: complex, gap_um: float = GAP_UM) -> float:
    """Потери на разрыве для заданной ширины полоски у края разрыва."""
    x = np.linspace(-40.0, 40.0, 1024)
    z = np.linspace(-40.0, 40.0, 1024)
    field, _ = mode_field(width_at_gap_um, x, z, n_planar)
    propagated = propagate_angular_spectrum(field, x, z, gap_um, N_CLAD, LAMBDA_UM)
    eta = overlap_2d(propagated, field, x, z)
    return coupling_loss_db(eta)


def gap_length_sweep(n_planar: complex) -> list[dict]:
    rows = []
    for g in (1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0):
        rows.append({"gap": g, "loss_db": gap_loss(W_GUIDE_UM, n_planar, g)})
    return rows


def taper_sweep(n_planar: complex) -> list[dict]:
    """Потери разрыва в зависимости от ширины полоски у его края."""
    rows = []
    for w in TAPER_WIDTHS:
        rows.append({
            "width": float(w),
            "gap_db": gap_loss(float(w), n_planar),
            "measured_db": MEASURED.get(float(w)),
        })
    return rows


def mode_width_from_gap_loss(target_db: float, gap_um: float = GAP_UM) -> float:
    """Обратная задача: какой размер моды даёт заданные потери на разрыве.

    Мода приближается гауссовым пучком, поскольку из измерения извлекается одна
    величина - её эффективный поперечный размер.
    """
    from slabmodes.coupling import gaussian_field
    x = np.linspace(-40.0, 40.0, 1024)
    z = np.linspace(-40.0, 40.0, 1024)

    def loss(mfd: float) -> float:
        f = gaussian_field(x, z, mfd)
        p = propagate_angular_spectrum(f, x, z, gap_um, N_CLAD, LAMBDA_UM)
        return coupling_loss_db(overlap_2d(p, f, x, z))

    lo, hi = 1.0, 30.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if loss(mid) > target_db:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def mode_width_from_propagation_loss(n_planar: complex, measured_db_cm: float,
                                     width_um: float = W_GUIDE_UM) -> float:
    """Оценка ширины моды по измеренным погонным потерям.

    Планарная мода бесконечной ширины лежит целиком над металлом. У полоски
    потери меньше во столько раз, какая доля мощности приходится на металл,
    поэтому доля = измеренные потери / планарные, а поперечный размер моды
    оценивается как ширина полоски, делённая на эту долю.
    """
    planar = propagation_loss_db_per_cm(n_planar, LAMBDA_UM)
    fraction = measured_db_cm / planar
    return float(width_um / fraction) if fraction > 0 else float("nan")


def make_plots(taper_rows: list[dict], gap_rows: list[dict], n_planar: complex) -> None:
    fig, axs = plt.subplots(1, 3, figsize=(15.6, 4.9))

    ax = axs[0]
    w = [r["width"] for r in taper_rows]
    ax.plot(w, [r["gap_db"] for r in taper_rows], "o-", ms=5, color="#0B6E99",
            label="расчёт: дифракция в разрыве")
    meas = [(r["width"], r["measured_db"]) for r in taper_rows if r["measured_db"]]
    ax.plot([m[0] for m in meas], [m[1] for m in meas], "k*", ms=14, label="измерено")
    ax.set_xlabel("ширина полоски у края разрыва, мкм")
    ax.set_ylabel("потери, дБ")
    ax.set_title(f"Разрыв {GAP_UM:g} мкм: роль обратного тейпера")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, frameon=False)

    ax = axs[1]
    g = [r["gap"] for r in gap_rows]
    ax.plot(g, [r["loss_db"] for r in gap_rows], "o-", ms=5, color="#B23A48")
    ax.axvline(GAP_UM, color="0.4", ls="--", lw=1.2, label=f"разрыв статьи {GAP_UM:g} мкм")
    ax.set_xlabel("длина разрыва, мкм")
    ax.set_ylabel("потери, дБ")
    ax.set_title(f"Потери против длины разрыва, полоска {W_GUIDE_UM:g} мкм")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, frameon=False)

    ax = axs[2]
    x = np.linspace(-25.0, 25.0, 1024)
    z = np.linspace(-25.0, 25.0, 1024)
    for width, style in ((2.0, "-"), (6.0, "--")):
        field, _ = mode_field(width, x, z, n_planar)
        prop = propagate_angular_spectrum(field, x, z, GAP_UM, N_CLAD, LAMBDA_UM)
        mid = len(z) // 2
        ax.plot(x, np.abs(field[:, mid]) / np.max(np.abs(field[:, mid])), style,
                lw=1.6, color="#0B6E99", label=f"до разрыва, {width:g} мкм")
        ax.plot(x, np.abs(prop[:, mid]) / np.max(np.abs(field[:, mid])), style,
                lw=1.6, color="#B23A48", label=f"после {GAP_UM:g} мкм, {width:g} мкм")
    ax.set_xlim(-25, 25)
    ax.set_xlabel("поперечная координата, мкм")
    ax.set_ylabel("|поле|, норм.")
    ax.set_title("Расплывание пучка в разрыве")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7.5, frameon=False)

    fig.suptitle("Lee & Lee 2019: туннелирование LR-SPP через разрыв волновода", fontsize=12)
    fig.tight_layout()
    fig.savefig(OUT / "lee2019_gap.png", dpi=170)
    plt.close(fig)


def main() -> int:
    n_planar = solve_mode(stack(), K0, complex(1.4512, 3e-5))
    n_strip = solve_strip(K0, n_planar**2, EPS_CLAD, W_GUIDE_UM)

    taper_rows = taper_sweep(n_planar)
    gap_rows = gap_length_sweep(n_planar)
    make_plots(taper_rows, gap_rows, n_planar)

    lines: list[str] = []
    add = lines.append
    add("Lee & Lee, J. Nanosci. Nanotechnol. 19, 6192 (2019) - аналитическое воспроизведение")
    add(f"Полоска золота {T_AU_UM * 1000:.0f} нм x {W_GUIDE_UM:g} мкм в полимере ZPU450, "
        f"разрыв {GAP_UM:g} мкм")
    add("")
    add("1. Мода волновода")
    add(f"   планарный предел  n_eff = {n_planar.real:.6f}   "
        f"потери {propagation_loss_db_per_cm(n_planar, LAMBDA_UM):7.2f} дБ/см")
    add(f"   полоска {W_GUIDE_UM:g} мкм    n_eff = {n_strip.real:.6f}   "
        f"потери {propagation_loss_db_per_cm(n_strip, LAMBDA_UM):7.2f} дБ/см"
        f"   (измерено {PAPER_PROP_DB_CM:g})")
    add("")
    add("2. Потери разрыва против ширины полоски у его края")
    add("   W, мкм   расчёт, дБ   измерено, дБ   отличие")
    for r in taper_rows:
        m = r["measured_db"]
        diff = f"{r['gap_db'] - m:+7.3f}" if m else "      -"
        add(f"   {r['width']:5.1f}   {r['gap_db']:10.3f}   "
            f"{(f'{m:12.3f}' if m else '           -')}   {diff}")
    base = next(r for r in taper_rows if r["width"] == 2.0)
    best = min(taper_rows, key=lambda r: r["gap_db"])
    add(f"   выигрыш от расширения 2 -> {best['width']:g} мкм: "
        f"{base['gap_db'] - best['gap_db']:.3f} дБ (измерено 0.27-0.29)")
    add("")
    add("3. Потери против длины разрыва (полоска 2 мкм)")
    add("   g, мкм   расчёт, дБ")
    for r in gap_rows:
        add(f"   {r['gap']:5.1f}   {r['loss_db']:10.3f}")
    add("   статья, расчёт предыдущей работы группы: около 0.03 дБ при разрыве 1 мкм")
    add("")
    add("4. Обратная задача: какой размер моды отвечает измерениям")
    mfd_gap = mode_width_from_gap_loss(MEASURED[2.0])
    mfd_taper = mode_width_from_gap_loss(MEASURED[6.0])
    mfd_loss = mode_width_from_propagation_loss(n_planar, PAPER_PROP_DB_CM)
    eim_mfd = 2.0 / (propagation_loss_db_per_cm(n_strip, LAMBDA_UM)
                     / propagation_loss_db_per_cm(n_planar, LAMBDA_UM))
    add(f"   из измеренных потерь разрыва без тейперов ({MEASURED[2.0]:.2f} дБ): "
        f"размер моды {mfd_gap:.2f} мкм")
    add(f"   из измеренных потерь разрыва с тейперами ({MEASURED[6.0]:.2f} дБ):  "
        f"размер моды {mfd_taper:.2f} мкм")
    add(f"   из измеренных погонных потерь ({PAPER_PROP_DB_CM:g} дБ/см):          "
        f"размер моды {mfd_loss:.2f} мкм")
    add(f"   тот же способ, применённый к расчёту ЭДП ({propagation_loss_db_per_cm(n_strip, LAMBDA_UM):.2f} дБ/см): "
        f"{eim_mfd:.1f} мкм")
    add("")
    add("   Обе измеренные величины формально указывают на моду около 3 мкм, тогда как")
    add("   ЭДП даёт примерно втрое более широкую. Однако полновекторный расчёт")
    add("   методом конечных элементов (scripts/fem_check_femwell.py) ЭДП здесь")
    add("   ПОДТВЕРЖДАЕТ: для полоски 6 мкм он даёт длину спадания поля 10 мкм, то есть")
    add("   мода действительно широкая. Значит, узкую моду из измерений вывести нельзя,")
    add("   и обратная задача выше говорит не о размере моды, а о том, что измеренные")
    add("   величины содержат не только тот механизм, который в них подставлен.")
    add("")
    add("   Отсюда основной вывод по этой статье: дифракция в разрыве даёт лишь")
    add(f"   {taper_rows[0]['gap_db']:.2f} дБ из измеренных {MEASURED[2.0]:.2f} дБ, то есть около десятой части.")
    add("   Оставшееся приходится на рассеяние на двух резких окончаниях металла,")
    add("   ограничивающих разрыв, а не на расплывание пучка. Это согласуется и с")
    add("   формулировкой самих авторов про преобразование размера моды с малым")
    add("   излучением: выигрыш обратных тейперов - это смягчение окончания металла,")
    add("   а не уменьшение дифракции.")
    add(f"   Для сведения: если бы всё измеренное было дифракцией, выигрыш тейперов")
    add(f"   отвечал бы расширению моды с {mfd_gap:.2f} до {mfd_taper:.2f} мкм, то есть")
    add(f"   всего на {100 * (mfd_taper / mfd_gap - 1):.0f} %, что для расширения полоски с 2 до 6 мкм неправдоподобно мало")
    add("   и служит дополнительным доводом против дифракционного объяснения.")

    with (OUT / "lee2019_gap.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["taper_width_um", "gap_loss_db_calc", "gap_loss_db_meas"])
        for r in taper_rows:
            w.writerow([f"{r['width']:.1f}", f"{r['gap_db']:.4f}",
                        f"{r['measured_db']:.3f}" if r["measured_db"] else ""])
        w.writerow([])
        w.writerow(["gap_length_um", "loss_db_calc"])
        for r in gap_rows:
            w.writerow([f"{r['gap']:.1f}", f"{r['loss_db']:.4f}"])

    (OUT / "lee2019_summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    for line in lines:
        print(line)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
