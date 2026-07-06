"""
Sensitivity sweep по нелинейной добавке Delta_eps_NL для Au наночастиц
в композите (ветка A из Obsidian-заметки "Нелинейность наночастиц Au для
EMA и расчетов композита", раздел "Формульный блок для процедуры Лерера").

Физическая модель (феноменологическая chi^(3)-ветка):
    eps_p(I, lam) = eps_L(lam) + Delta_eps_NL,
    Delta_eps_NL  = (3/2) * chi3 * |f|^2 * I / (n_h * eps0 * c),
    f             = 3*eps_h / (eps_p + 2*eps_h)   (квазистатический фактор
                                                   локального поля сферы),
далее eps_eff(I, lam) через composite_from_eps() (формула КМ/МГ, как в
COMPOSIT.c). Расчет самосогласованный: f зависит от eps_p, eps_p - от
|E_loc|^2; решается итерациями с демпфированием.

Запись Delta_eps через интенсивность конвенционно-инвариантна: она одинакова
для амплитудных конвенций E(t)=(1/2)(E e^{-iwt}+к.с.) (Delta_eps=(3/4)chi|E|^2)
и E(t)=E e^{-iwt}+к.с. (Delta_eps=3chi|E|^2, как у De Leon 2014). Для
esu-значений остается фактор-4 неопределенность конвенции первоисточника -
на фоне межлитературного разброса chi^(3) в порядки она несущественна для
sweep (см. заметку).

Конвенция знака времени: таблицы этого порта хранят eps = eps' - i*eps''
(конвенция Лерера, e^{+i w t}, пассивная среда Im eps <= 0), поэтому
литературные chi^(3) (оптическая конвенция e^{-i w t}, Im eps >= 0)
КОМПЛЕКСНО СОПРЯГАЮТСЯ перед подстановкой: chi_Lerer = conj(chi_opt).

Статус: sensitivity study / diagnostic. Значения chi^(3) - literature
benchmarks, сверенные по PDF (см. Obsidian-заметку, раздел "Сверка
benchmarks по PDF (2026-07-06)"); они относятся к разным режимам
(нс/пс/фс, разные lambda и матрицы) и НЕ являются откалиброванными
параметрами конкретного образца. Цель - определить, при каких I и chi^(3)
эффект в eps_eff становится заметным, а где ломается возмущательная ветка.

Benchmarks (оптическая конвенция, до сопряжения):
- Hache 1988 (Appl. Phys. A 47, 347): |chi_m| ~ 5e-8 esu, фаза ~80 град
  (преимущественно мнимая, hot-electron), 527/532 нм, 5 ps,
  I <= 4-5e8 Вт/см^2. PDF сверен.
- Guillet 2009 (PRB 79, 045410, теория/TTM): |chi_he| = 2.9e-7 esu,
  фаза 84 град, hbar*w=1.5 эВ, 2 ps. PDF сверен.
- De Leon 2014 (Opt. Lett. 39, 2274): chi_Au = (4.67+3.03i)e-19 м^2/В^2,
  795.6 нм, 107 фс, пленка Au (SPP). PDF сверен.
- Rotenberg (по пересчету De Leon 2014): chi ~ (-7.6+0.4i)e-19 м^2/В^2,
  630 нм, 100 фс - знак Re меняется с lambda. Вторичный пересчет.
- Ryasnyanskiy 2007 (J. Lumin. 127, 181): |chi_m| = 6.25e-8 esu (Au:Al2O3,
  532 нм, 7 нс) - в статье только МОДУЛЬ интринсик-восприимчивости, фаза
  неизвестна; здесь брекетинг фаз 0/90/180/270 град. PDF сверен.

Пределы применимости по I (из тех же PDF): нс-эксперименты 5.7-28 МВт/см^2;
пс - до ~5e8 Вт/см^2; фс Z-scan 2e10-1.4e11 Вт/см^2; насыщение SA сфер
~6e10 Вт/см^2 (Ahmedov 2026); оптический пробой ~4e12 Вт/см^2 (200 фс).

Запуск (из корня репозитория):
    python composite_ema/scripts/nonlinear_au_sweep.py
Выход: composite_ema/results/nonlinear_au_sweep_{spectra,intensity}.{csv,png},
       nonlinear_au_sweep_summary.csv + печать самопроверок.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from composite_ema import (
    RESULTS_DIR,
    composite_from_eps,
    eps_au,
    find_frohlich_crossing,
)

EPS0 = 8.8541878128e-12  # Ф/м
C_LIGHT = 2.99792458e8   # м/с
ESU_TO_SI = 1.3963e-8    # chi^(3): 1 esu -> м^2/В^2 (Boyd: 4*pi/(3e4)^2)

N_HOST = 1.77            # показатель преломления слоя с ПНЧ (как в статье)
C_FILL = 0.10            # объемная концентрация ПНЧ (как в статье)


@dataclass(frozen=True)
class ChiBenchmark:
    """Литературный benchmark chi^(3) в ОПТИЧЕСКОЙ конвенции e^{-i w t}."""

    name: str
    chi_opt_si: complex     # м^2/В^2, оптическая конвенция
    lam_ref_nm: float       # длина волны первоисточника
    i_max_wcm2: float       # максимальная интенсивность режима первоисточника
    note: str

    @property
    def chi_lerer_si(self) -> complex:
        """chi^(3) в конвенции Лерера e^{+i w t}: комплексное сопряжение."""
        return self.chi_opt_si.conjugate()


def _from_esu(modulus_esu: float, phase_deg: float) -> complex:
    phase = math.radians(phase_deg)
    return modulus_esu * ESU_TO_SI * complex(math.cos(phase), math.sin(phase))


BENCHMARKS = [
    ChiBenchmark(
        "Hache1988_ps_SPR",
        _from_esu(5.0e-8, 80.0),
        527.0, 5.0e8,
        "Au сферы в стекле, OPC/DFWM 5 ps; |chi_m|~5e-8 esu, фаза ~80 град (hot-electron)",
    ),
    ChiBenchmark(
        "Guillet2009_2ps_TTM",
        _from_esu(2.9e-7, 84.0),
        827.0, 5.0e8,
        "теория TTM, hbar*w=1.5 эВ, 2 ps; верхняя пс-оценка hot-electron ветки",
    ),
    ChiBenchmark(
        "DeLeon2014_fs_795nm",
        complex(4.67e-19, 3.03e-19),
        795.6, 1.0e11,
        "пленка Au (SPP, Кречман), 107 фс; прямое комплексное chi_Au",
    ),
    ChiBenchmark(
        "Rotenberg_fs_630nm",
        complex(-7.6e-19, 0.4e-19),
        630.0, 1.0e11,
        "пересчет De Leon 2014 из Z-scan Rotenberg (n2=0); знак Re отрицательный",
    ),
]

# Ryasnyanskiy 2007: известен только модуль intrinsic chi_m - брекетинг фаз.
RYASN_MODULUS_ESU = 6.25e-8
RYASN_PHASES_DEG = (0.0, 90.0, 180.0, 270.0)
BENCHMARKS_RYASN = [
    ChiBenchmark(
        f"Ryasn2007_ns_532_phase{int(ph)}",
        _from_esu(RYASN_MODULUS_ESU, ph),
        532.0, 2.8e7,
        f"Au:Al2O3 7 ns; |chi_m|=6.25e-8 esu, фаза {int(ph)} град (брекетинг, фаза в статье не измерена)",
    )
    for ph in RYASN_PHASES_DEG
]

# Отметки применимости по интенсивности (Вт/см^2), сверены по PDF:
I_MARKS_WCM2 = {
    "ns Z-scan (Ryasnyanskiy)": 2.8e7,
    "ps OPC (Hache)": 5.0e8,
    "SA сфер, Isat (Ahmedov)": 6.0e10,
    "пробой 200 fs (Ahmedov)": 4.0e12,
}


def local_field_factor(eps_p: complex, eps_h: complex) -> complex:
    return 3.0 * eps_h / (eps_p + 2.0 * eps_h)


def eps_p_nonlinear(
    lam_nm: float,
    i_wcm2: float,
    chi_lerer_si: complex,
    eps_h: complex,
    max_iter: int = 200,
    tol: float = 1e-12,
) -> tuple[complex, bool]:
    """Самосогласованное eps_p(I): Delta_eps = (3/2)chi*|f|^2*I/(n_h eps0 c).

    Возвращает (eps_p, converged). Демпфирование 0.5 для устойчивости у
    резонанса. I подается в Вт/см^2 и переводится в Вт/м^2 внутри.
    """
    eps_l = eps_au(lam_nm)
    if i_wcm2 == 0.0:
        return eps_l, True
    i_si = i_wcm2 * 1.0e4
    n_h = math.sqrt(eps_h.real)
    prefactor = 1.5 * chi_lerer_si * i_si / (n_h * EPS0 * C_LIGHT)
    eps_p = eps_l
    for _ in range(max_iter):
        f = local_field_factor(eps_p, eps_h)
        eps_new = eps_l + prefactor * (abs(f) ** 2)
        if abs(eps_new - eps_p) < tol * max(1.0, abs(eps_new)):
            return eps_new, True
        eps_p = 0.5 * eps_p + 0.5 * eps_new
    return eps_p, False


def perturbation_status(eps_p: complex, lam_nm: float) -> str:
    """Классификация достоверности линейной по I ветки."""
    delta = abs(eps_p - eps_au(lam_nm))
    scale = abs(eps_au(lam_nm))
    if delta <= 0.1 * scale:
        return "ok"
    if delta <= 1.0 * scale:
        return "questionable"
    return "breakdown"


def self_test() -> None:
    eps_h = complex(N_HOST**2, 0.0)

    # 1) I=0 возвращает ровно линейный композит.
    for lam in (450.0, 532.0, 650.0, 795.6):
        eps_p0, conv = eps_p_nonlinear(lam, 0.0, BENCHMARKS[0].chi_lerer_si, eps_h)
        assert conv and eps_p0 == eps_au(lam)
        assert composite_from_eps(eps_h, eps_p0, C_FILL) == composite_from_eps(
            eps_h, eps_au(lam), C_FILL
        )
    print("[self_test] OK: I=0 воспроизводит линейный eps_eff без изменений.")

    # 2) Дилютный аналитический контроль: d(eps_eff) = C f^2 Delta_eps_p
    #    (в пределе C->0), т.е. chi_eff = C f^2 |f|^2 chi_m; точная производная
    #    формулы МГ дает дополнительный множитель 1/(1-C*eta)^2.
    lam = 532.0
    chi = BENCHMARKS[0].chi_lerer_si
    i_small = 1.0e3  # Вт/см^2 - глубоко линейный режим
    # допуск дилютной формулы ~ |2*C*eta|: у резонанса |eta|~2.9, поэтому
    # уже при C=0.01 отклонение ~6%; точная поправка (1-C*eta)^-2 обязана
    # сходиться до <0.1% при любом C (второй assert).
    for c_fill, tol in ((0.01, 0.10), (C_FILL, 0.75)):
        eps_l = eps_au(lam)
        eps_p, conv = eps_p_nonlinear(lam, i_small, chi, eps_h)
        assert conv
        d_eff_num = composite_from_eps(eps_h, eps_p, c_fill) - composite_from_eps(
            eps_h, eps_l, c_fill
        )
        f_lin = local_field_factor(eps_l, eps_h)
        d_eff_ana = c_fill * f_lin**2 * (eps_p - eps_l)
        rel = abs(d_eff_num - d_eff_ana) / abs(d_eff_num)
        # точная поправка (1 - C*eta)^{-2}:
        eta = (eps_l - eps_h) / (eps_l + 2.0 * eps_h)
        d_eff_exact = d_eff_ana / (1.0 - c_fill * eta) ** 2
        rel_exact = abs(d_eff_num - d_eff_exact) / abs(d_eff_num)
        assert rel < tol, (c_fill, rel)
        assert rel_exact < 1e-3, (c_fill, rel_exact)
        print(
            f"[self_test] OK: дилютный контроль C={c_fill:.2f}: "
            f"chi_eff=C f^2|f|^2 chi_m с отклонением {rel:.1%} "
            f"(с поправкой (1-C*eta)^-2: {rel_exact:.2%})."
        )

    # 3) Кросс-чек Ryasnyanskiy 2007: |chi_eff| = p |f|^4 |chi_m| на 532 нм
    #    против измеренного composite-уровня |chi^(3)| = 7.33e-7 esu
    #    (Au:Al2O3, p=8%, n_host Al2O3 ~ 1.77). Наши eps Au - J&C 1972,
    #    у авторов другая таблица, поэтому ждем совпадения по порядку.
    f532 = local_field_factor(eps_au(532.0), eps_h)
    chi_eff_esu = 0.08 * abs(f532) ** 4 * RYASN_MODULUS_ESU
    ratio = chi_eff_esu / 7.33e-7
    assert 0.2 < ratio < 5.0, ratio
    print(
        f"[self_test] OK: |chi_eff| по формуле p|f|^4|chi_m| = {chi_eff_esu:.2e} esu; "
        f"измерено Ryasnyanskiy 7.33e-7 esu (отношение {ratio:.2f} - порядок совпадает)."
    )

    # 4) Знак: сопряженная hot-electron ветка (Im chi_opt > 0) в конвенции
    #    Лерера должна ДОБАВЛЯТЬ потери частице: Im eps_p уходит вниз
    #    (более отрицательный), пассивность частицы сохраняется.
    eps_p, conv = eps_p_nonlinear(532.0, 1.0e6, BENCHMARKS[0].chi_lerer_si, eps_h)
    assert conv and eps_p.imag < eps_au(532.0).imag <= 0.0
    print(
        "[self_test] OK: hot-electron benchmark в конвенции Лерера увеличивает "
        "|Im eps_p| (потери частицы растут, Im eps_p <= 0 сохраняется)."
    )


def run_spectral_sweep() -> None:
    """Спектры eps_eff(lambda) при нескольких I для Hache-benchmark."""
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    eps_h = complex(N_HOST**2, 0.0)
    bench = BENCHMARKS[0]  # Hache 1988 - пс-режим на SPR, ближе всего к КМ-слою
    intensities = [0.0, 1.0e6, 5.0e6, 2.0e7]  # Вт/см^2, зона нс/пс экспериментов
    lam_grid = np.arange(400.0, 900.0 + 1e-9, 2.0)

    rows = []
    curves: dict[float, tuple[list[float], list[float], list[float]]] = {
        i: ([], [], []) for i in intensities
    }
    for lam in lam_grid:
        row = {"lambda_nm": lam}
        for i_wcm2 in intensities:
            eps_p, conv = eps_p_nonlinear(lam, i_wcm2, bench.chi_lerer_si, eps_h)
            eps_eff = composite_from_eps(eps_h, eps_p, C_FILL)
            status = perturbation_status(eps_p, lam)
            key = f"I{i_wcm2:.0e}"
            row[f"Re_eps_eff_{key}"] = eps_eff.real
            row[f"Im_eps_eff_{key}"] = eps_eff.imag
            row[f"status_{key}"] = status if conv else "no_convergence"
            lam_l, re_l, im_l = curves[i_wcm2]
            lam_l.append(lam)
            re_l.append(eps_eff.real)
            im_l.append(eps_eff.imag)
        rows.append(row)

    csv_path = RESULTS_DIR / "nonlinear_au_sweep_spectra.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"[spectra] CSV записан: {csv_path}")

    cross = find_frohlich_crossing(eps_au, eps_h.real, 400.0, 900.0)
    print(f"[spectra] Условие Фрёлиха Au (n_host={N_HOST}): {cross.lam_nm:.1f} нм")

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax_re, ax_im) = plt.subplots(2, 1, figsize=(7.5, 7.5), sharex=True)
    colors = plt.cm.viridis(np.linspace(0.0, 0.85, len(intensities)))
    for color, i_wcm2 in zip(colors, intensities):
        lam_l, re_l, im_l = curves[i_wcm2]
        label = "линейный (I=0)" if i_wcm2 == 0 else f"I = {i_wcm2:.0e} Вт/см²"
        ax_re.plot(lam_l, re_l, color=color, label=label)
        ax_im.plot(lam_l, im_l, color=color, label=label)
    if cross.lam_nm:
        for ax in (ax_re, ax_im):
            ax.axvline(cross.lam_nm, color="gray", linestyle=":", linewidth=1)
    ax_re.set_ylabel(r"$\mathrm{Re}(\varepsilon_\mathrm{eff})$")
    ax_re.legend(fontsize=8)
    ax_re.set_title(
        f"Delta_eps_NL sweep (Hache 1988, 5e-8 esu @80°): n_host={N_HOST}, C={C_FILL:.0%}",
        fontsize=10,
    )
    ax_im.set_xlabel("длина волны, нм")
    ax_im.set_ylabel(r"$\mathrm{Im}(\varepsilon_\mathrm{eff})$")
    ax_im.legend(fontsize=8)
    fig.tight_layout()
    png_path = RESULTS_DIR / "nonlinear_au_sweep_spectra.png"
    fig.savefig(png_path, dpi=150)
    print(f"[spectra] PNG записан: {png_path}")


def run_intensity_sweep() -> None:
    """delta(I) = |Im eps_eff(I) - Im eps_eff(0)| / |Im eps_eff(0)| на резонансе."""
    eps_h = complex(N_HOST**2, 0.0)
    cross = find_frohlich_crossing(eps_au, eps_h.real, 400.0, 900.0)
    lam_res = cross.lam_nm if cross.lam_nm else 550.0
    i_grid = np.logspace(4.0, 12.0, 49)  # Вт/см^2

    all_bench = BENCHMARKS + BENCHMARKS_RYASN
    eps_eff_lin = composite_from_eps(eps_h, eps_au(lam_res), C_FILL)

    rows = []
    summary = []
    curves: dict[str, tuple[list[float], list[float], list[str]]] = {}
    for bench in all_bench:
        i_list: list[float] = []
        d_list: list[float] = []
        s_list: list[str] = []
        i_at_1pct = None
        i_pert_limit = None
        for i_wcm2 in i_grid:
            eps_p, conv = eps_p_nonlinear(lam_res, float(i_wcm2), bench.chi_lerer_si, eps_h)
            status = perturbation_status(eps_p, lam_res) if conv else "no_convergence"
            eps_eff = composite_from_eps(eps_h, eps_p, C_FILL)
            delta_rel = abs(eps_eff.imag - eps_eff_lin.imag) / abs(eps_eff_lin.imag)
            if eps_eff.imag > 0.0:
                status = "breakdown"  # композит перестал быть пассивным
            rows.append(
                {
                    "benchmark": bench.name,
                    "lambda_nm": lam_res,
                    "I_wcm2": float(i_wcm2),
                    "Re_eps_eff": eps_eff.real,
                    "Im_eps_eff": eps_eff.imag,
                    "delta_Im_rel": delta_rel,
                    "status": status,
                }
            )
            i_list.append(float(i_wcm2))
            d_list.append(delta_rel)
            s_list.append(status)
            if i_at_1pct is None and delta_rel >= 0.01 and status == "ok":
                i_at_1pct = float(i_wcm2)
            if i_pert_limit is None and status != "ok":
                i_pert_limit = float(i_wcm2)
        curves[bench.name] = (i_list, d_list, s_list)
        summary.append(
            {
                "benchmark": bench.name,
                "chi_opt_si_re": bench.chi_opt_si.real,
                "chi_opt_si_im": bench.chi_opt_si.imag,
                "lam_ref_nm": bench.lam_ref_nm,
                "I_source_max_wcm2": bench.i_max_wcm2,
                "I_at_1pct_dImEff_wcm2": i_at_1pct,
                "I_perturbation_limit_wcm2": i_pert_limit,
                "note": bench.note,
            }
        )

    csv_path = RESULTS_DIR / "nonlinear_au_sweep_intensity.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"[intensity] CSV записан: {csv_path}")

    sum_path = RESULTS_DIR / "nonlinear_au_sweep_summary.csv"
    with sum_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(summary[0].keys()))
        writer.writeheader()
        writer.writerows(summary)
    print(f"[intensity] Сводка записана: {sum_path}")
    print(f"[intensity] Резонансная длина волны: {lam_res:.1f} нм; "
          f"Im eps_eff(лин) = {eps_eff_lin.imag:.4f}")
    for s in summary:
        print(
            f"[intensity]   {s['benchmark']}: 1% изменения Im eps_eff при "
            f"I ~ {s['I_at_1pct_dImEff_wcm2']!s} Вт/см²; предел возмущательной "
            f"ветки I ~ {s['I_perturbation_limit_wcm2']!s} Вт/см² "
            f"(режим источника до {s['I_source_max_wcm2']:.0e})"
        )

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8.5, 6.0))
    for bench in all_bench:
        i_list, d_list, s_list = curves[bench.name]
        i_ok = [i for i, s in zip(i_list, s_list) if s == "ok"]
        d_ok = [d for d, s in zip(d_list, s_list) if s == "ok"]
        i_bad = [i for i, s in zip(i_list, s_list) if s != "ok"]
        d_bad = [d for d, s in zip(d_list, s_list) if s != "ok"]
        style = "--" if bench.name.startswith("Ryasn") else "-"
        (line,) = ax.loglog(i_ok, d_ok, style, label=bench.name, linewidth=1.4)
        if i_bad:
            ax.loglog(i_bad, d_bad, ":", color=line.get_color(), alpha=0.45, linewidth=1.0)
    ax.axhline(0.01, color="black", linewidth=0.8, linestyle="-.")
    ax.text(1.3e4, 0.011, "1% изменения Im eps_eff", fontsize=8)
    for label, i_mark in I_MARKS_WCM2.items():
        ax.axvline(i_mark, color="gray", linestyle=":", linewidth=0.8)
        ax.text(i_mark * 1.1, 2e-7, label, rotation=90, fontsize=7, color="gray")
    ax.set_xlabel("I, Вт/см²")
    ax.set_ylabel(r"$|\Delta \mathrm{Im}\,\varepsilon_\mathrm{eff}| / |\mathrm{Im}\,\varepsilon_\mathrm{eff,lin}|$")
    ax.set_title(
        f"Чувствительность eps_eff к Delta_eps_NL на {lam_res:.0f} нм "
        f"(n_host={N_HOST}, C={C_FILL:.0%});\nпунктирное продолжение = вне возмущательной ветки",
        fontsize=9,
    )
    ax.set_ylim(1e-7, 30)
    ax.legend(fontsize=7, loc="upper left")
    ax.grid(True, which="both", alpha=0.25)
    fig.tight_layout()
    png_path = RESULTS_DIR / "nonlinear_au_sweep_intensity.png"
    fig.savefig(png_path, dpi=150)
    print(f"[intensity] PNG записан: {png_path}")


if __name__ == "__main__":
    self_test()
    run_spectral_sweep()
    run_intensity_sweep()
