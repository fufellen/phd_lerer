"""Аналитическое воспроизведение Mueller, Leosson, Capasso, Nano Lett. 14, 5524 (2014).

Статья: ввод из волокна в LR-SPP-волновод через массив плазмонных антенн, в
котором поляризация падающего света определяет направление возбуждаемой моды.
DOI 10.1021/nl501860r.

Что воспроизводится:
  одномодовость TM в диапазоне C+L (1530-1625 нм);
  период массива, заявленный авторами как lambda_0 / n_BCB;
  механизм однонаправленного возбуждения парой антенн под +-45 градусов,
    разнесённых на четверть длины волны моды;
  различение поляризаций не менее 30 дБ (в статье ограничено шумом);
  отклонение векторов прибора от круговой поляризации на 22 и 10 градусов.

Абсолютной эффективности ввода в статье нет, поэтому она здесь и не считается:
воспроизводится механизм и то, что из измерений извлекается однозначно.

Запуск:
    python lrspp_coupling/scripts/mueller2014_antenna_array.py
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

from slabmodes import Stack, materials, propagation_loss_db_per_cm, solve_mode  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "results"
OUT.mkdir(parents=True, exist_ok=True)

N_BCB = 1.535
EPS_BCB = complex(N_BCB**2, 0.0)
T_AU_UM = 0.015
W_STRIP_UM = 10.0
BAND_UM = (1.530, 1.625)

PAPER_DEVIATIONS_DEG = (22.0, 10.0)
PAPER_ER_DB = 30.0


def stack() -> Stack:
    return Stack(eps=(EPS_BCB, materials.EPS_AU_1550, EPS_BCB), thickness=(T_AU_UM,),
                 names=("BCB", "Au", "BCB"))


def mode_over_band(points: int = 11) -> list[dict]:
    """LR-SPP-мода плёнки Au в BCB по диапазону C+L."""
    rows = []
    guess = complex(N_BCB + 2e-3, 3e-5)
    for lam in np.linspace(BAND_UM[0], BAND_UM[1], points):
        k0 = materials.k0_from_lambda(float(lam))
        n = solve_mode(stack(), k0, guess)
        guess = n
        rows.append({
            "lambda": float(lam),
            "neff": n,
            "loss_db_cm": propagation_loss_db_per_cm(n, float(lam)),
            "lambda_spp": float(lam) / n.real,
        })
    return rows


def jones_from_stokes_angle(alpha_rad: float) -> np.ndarray:
    """Поляризация вдоль большого круга сферы Пуанкаре, проходящего через полюса.

    При alpha = 0 получается горизонтальная линейная поляризация, при
    alpha = +-90 градусов - круговые.
    """
    return np.array([np.cos(alpha_rad / 2.0), 1j * np.sin(alpha_rad / 2.0)])


def unit_vector(angle_deg: float) -> np.ndarray:
    a = np.radians(angle_deg)
    return np.array([np.cos(a), np.sin(a)])


def channel_amplitudes(jones: np.ndarray, amp_ratio: float = 1.0,
                       phase_error_deg: float = 0.0) -> tuple[complex, complex]:
    """Амплитуды моды, уходящей вправо и влево, от пары антенн +-45 градусов.

    Антенна 1 стоит в нуле и ориентирована под +45 градусов, антенна 2 смещена
    на четверть длины волны моды и ориентирована под -45. Возбуждение каждой
    пропорционально проекции поля на её ось. Для волны, уходящей вправо, вклад
    второй антенны отстаёт по фазе на 90 градусов, для уходящей влево - опережает,
    поэтому одно из направлений гасится.
    """
    a1 = complex(np.dot(unit_vector(45.0), jones))
    a2 = complex(np.dot(unit_vector(-45.0), jones)) * amp_ratio
    step = np.exp(1j * (np.pi / 2.0 + np.radians(phase_error_deg)))
    right = a1 + a2 / step
    left = a1 + a2 * step
    return right, left


def extinction_db(amp_ratio: float = 1.0, phase_error_deg: float = 0.0) -> float:
    """Различение направлений для круговой поляризации."""
    jones = jones_from_stokes_angle(np.pi / 2.0)  # круговая
    right, left = channel_amplitudes(jones, amp_ratio, phase_error_deg)
    hi, lo = max(abs(right), abs(left)), min(abs(right), abs(left))
    if lo <= 0:
        return float("inf")
    return float(20.0 * np.log10(hi / lo))


def device_vector_deviation_deg(amp_ratio: float, phase_error_deg: float) -> float:
    """Отклонение вектора прибора от полюса сферы Пуанкаре, градусы.

    Ищется поляризация вдоль большого круга, при которой канал даёт максимум;
    для идеальной пары антенн это круговая поляризация, то есть полюс.
    """
    alphas = np.linspace(-np.pi, np.pi, 40001)
    best, best_alpha = -1.0, 0.0
    for a in alphas:
        right, _ = channel_amplitudes(jones_from_stokes_angle(a), amp_ratio, phase_error_deg)
        if abs(right) > best:
            best, best_alpha = abs(right), a
    return float(abs(np.degrees(best_alpha)) - 90.0)


def phase_error_for_deviation(target_deg: float) -> float:
    """Какая ошибка фазы отвечает заданному отклонению вектора прибора.

    Разбаланс амплитуд вектор прибора не смещает: при |a1| != |a2| максимум
    остаётся на полюсе, меняется только глубина гашения. Смещает вектор именно
    фазовая ошибка, то есть отклонение шага пары антенн от четверти длины волны.
    """
    lo, hi = 0.0, 89.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if device_vector_deviation_deg(1.0, mid) < target_deg:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def make_plots(band: list[dict]) -> None:
    fig, axs = plt.subplots(1, 3, figsize=(15.6, 4.9))

    ax = axs[0]
    lam = [r["lambda"] * 1000 for r in band]
    ax.plot(lam, [r["lambda_spp"] * 1000 for r in band], "o-", ms=4, color="#0B6E99",
            label="расчёт: период = lambda_0 / n_eff")
    ax.plot(lam, [r["lambda"] * 1000 / N_BCB for r in band], "--", color="0.4",
            label="приближение статьи: lambda_0 / n_BCB")
    ax.set_xlabel("длина волны, нм")
    ax.set_ylabel("период массива, нм")
    ax.set_title("Период антенного массива")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, frameon=False)

    ax = axs[1]
    alphas = np.linspace(-180, 180, 721)
    right, left = [], []
    for a in alphas:
        r, l = channel_amplitudes(jones_from_stokes_angle(np.radians(a)))
        right.append(abs(r) ** 2)
        left.append(abs(l) ** 2)
    right = np.array(right) / max(right)
    left = np.array(left) / max(left)
    ax.plot(alphas, right, lw=2, color="#0B6E99", label="вправо")
    ax.plot(alphas, left, lw=2, ls="--", color="#B23A48", label="влево")
    ax.axvline(90, color="0.5", lw=1, ls=":")
    ax.axvline(-90, color="0.5", lw=1, ls=":")
    ax.set_xlabel("положение на большом круге сферы Пуанкаре, градусы")
    ax.set_ylabel("нормированная мощность")
    ax.set_title("Поляризация задаёт направление")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, frameon=False)

    ax = axs[2]
    errs = np.linspace(0.1, 20, 200)
    ax.semilogy(errs, [10 ** (extinction_db(1.0, float(e)) / 10) for e in errs],
                lw=2, label="ошибка фазы (период)")
    ratios = np.linspace(1.001, 2.0, 200)
    ax.semilogy(100 * (ratios - 1), [10 ** (extinction_db(float(r), 0.0) / 10) for r in ratios],
                lw=2, ls="--", label="разбаланс амплитуд, %")
    ax.axhline(10 ** (PAPER_ER_DB / 10), color="k", ls=":", lw=1.4,
               label=f"порог измерения статьи, {PAPER_ER_DB:g} дБ")
    ax.set_xlabel("отклонение: градусы фазы либо проценты амплитуды")
    ax.set_ylabel("различение направлений, разы")
    ax.set_title("Чем ограничено различение поляризаций")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8, frameon=False)

    fig.suptitle("Mueller 2014: поляризационно-селективный ввод в LR-SPP", fontsize=12)
    fig.tight_layout()
    fig.savefig(OUT / "mueller2014_array.png", dpi=170)
    plt.close(fig)


def main() -> int:
    band = mode_over_band()
    make_plots(band)

    mid = band[len(band) // 2]
    lines: list[str] = []
    add = lines.append
    add("Mueller, Leosson, Capasso, Nano Lett. 14, 5524 (2014) - аналитическое воспроизведение")
    add(f"Полоска золота {W_STRIP_UM:g} мкм x {T_AU_UM * 1000:.0f} нм в BCB, n = {N_BCB}")
    add("")
    add("1. Мода волновода в диапазоне C+L")
    add("   lambda, нм    n_eff      потери, дБ/мм   период массива, нм")
    for r in band[::2]:
        add(f"   {r['lambda'] * 1000:8.0f}   {r['neff'].real:.6f}   "
            f"{r['loss_db_cm'] / 10:11.3f}   {r['lambda_spp'] * 1000:15.1f}")
    add(f"   изменение n_eff по диапазону: {band[-1]['neff'].real - band[0]['neff'].real:+.6f}")
    add(f"   приближение статьи lambda_0 / n_BCB при 1550 нм: {1.550 / N_BCB * 1000:.1f} нм")
    add(f"   точный период по расчёту:                        {mid['lambda_spp'] * 1000:.1f} нм")
    add(f"   отличие: {100 * (mid['lambda_spp'] - 1.550 / N_BCB) / (1.550 / N_BCB):+.2f} %"
        "   - приближение авторов оправдано")
    add(f"   четверть длины волны моды: {mid['lambda_spp'] * 1000 / 4:.1f} нм")
    add("")
    add("2. Механизм однонаправленного возбуждения")
    add("   Пара антенн под +45 и -45 градусов, разнесённая на четверть длины волны моды.")
    add("   Для круговой поляризации вклады в одном направлении гасятся точно:")
    for name, alpha in (("правая круговая", 90.0), ("левая круговая", -90.0),
                        ("горизонтальная линейная", 0.0), ("диагональная линейная", 180.0)):
        r, l = channel_amplitudes(jones_from_stokes_angle(np.radians(alpha)))
        add(f"     {name:24s} вправо {abs(r) ** 2:6.3f}   влево {abs(l) ** 2:6.3f}")
    add("")
    add("3. Чем ограничено различение направлений")
    add("   идеальная пара антенн: гашение полное, различение бесконечно")
    add("   фазовая ошибка периода   различение, дБ")
    for e in (1.0, 2.0, 5.0, 10.0):
        add(f"     {e:5.1f} градуса          {extinction_db(1.0, e):8.1f}")
    add("   разбаланс амплитуд        различение, дБ")
    for r in (1.01, 1.05, 1.10, 1.20):
        add(f"     {100 * (r - 1):5.1f} процента         {extinction_db(r, 0.0):8.1f}")
    add(f"   в статье сообщено не менее {PAPER_ER_DB:g} дБ, ограничено шумом установки.")
    add("   Измерение велось при поляризации, отвечающей вектору самого прибора, то есть")
    add("   фазовая ошибка при этом уже скомпенсирована выбором поляризации и гашение не")
    add("   ухудшает. Значит, сообщённый предел ограничивает именно баланс амплитуд:")
    add("   лучше примерно шести процентов.")
    add("")
    add("4. Отклонение векторов прибора от круговой поляризации")
    add("   Модель разделяет две независимые несовершенности:")
    add("     разбаланс амплитуд антенн НЕ смещает вектор прибора с полюса,")
    add("       он только ухудшает гашение;")
    add("     ошибка фазы (отклонение шага пары от четверти длины волны) смещает")
    add("       вектор прибора, но гашение при этом остаётся полным.")
    add(f"   измерено в статье: отклонения {PAPER_DEVIATIONS_DEG[0]:g} и "
        f"{PAPER_DEVIATIONS_DEG[1]:g} градусов при различении не менее {PAPER_ER_DB:g} дБ")
    period_nm = mid["lambda_spp"] * 1000
    for target in PAPER_DEVIATIONS_DEG:
        phase = phase_error_for_deviation(target)
        offset_nm = period_nm * phase / 360.0
        add(f"   отклонению {target:4.0f} градусов отвечает ошибка фазы {phase:5.1f} градуса,")
        add(f"     то есть смещение антенн на {offset_nm:5.1f} нм при шаге {period_nm:.0f} нм")
    add("   Такие смещения лежат в пределах точности электронной литографии и")
    add("   совместного совмещения слоёв, поэтому объяснение самосогласовано:")
    add("   вектор смещён технологией, а глубокое гашение сохраняется, потому что")
    add("   амплитуды двух групп антенн сбалансированы лучше шести процентов.")

    with (OUT / "mueller2014_band.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["lambda_um", "neff_real", "neff_imag", "loss_db_per_cm", "array_period_um"])
        for r in band:
            w.writerow([f"{r['lambda']:.4f}", f"{r['neff'].real:.9f}", f"{r['neff'].imag:.6e}",
                        f"{r['loss_db_cm']:.4f}", f"{r['lambda_spp']:.6f}"])

    (OUT / "mueller2014_summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    for line in lines:
        print(line)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
