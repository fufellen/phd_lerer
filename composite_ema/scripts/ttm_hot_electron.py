"""
Ветка B нелинейности Au ПНЧ: двухтемпературная модель (TTM) hot-electron /
thermal динамики. Диагностический расчет к задаче Лерера "нелинейность н/ч Au".

Модель (Guillet et al., PRB 79, 045410 (2009), уравнения (9)-(10)):
    C_e(T_e) dT_e/dt = -G (T_e - T_l) + P_abs(t),   C_e = gamma_e * T_e,
    C_l        dT_l/dt =  G (T_e - T_l),
где P_abs(t) - объемная плотность поглощенной мощности в Au (гауссов импульс).

Параметры золота (все ПЕРВИЧНЫЕ, из Guillet 2009, сверено по PDF 2026-07-06):
    gamma_e = 66      Дж·м^-3·К^-2   (C_e = gamma_e T_e; ссылка [22])
    C_l     = 2.44e6  Дж·м^-3·К^-1   (закон Дюлонга-Пти, T_l > T_D)
    T_D     = 170     К
    G       = 3e16    Вт·м^-3·К^-1   (electron-phonon coupling; ссылка [36])

Что делает скрипт (и чего НЕ делает):
- решает TTM для представительных импульсов, проверяет сохранение энергии,
  масштаб времени e-ph релаксации (~ps, сверка с Hua 2015: 1.8-5.7 ps) и
  переход к квазиравновесию T_e≈T_l при tau_p >> tau_ep (нс/CW-режим);
- НЕ пересчитывает Delta_eps из T_e, T_l: для этого нужны спектральные
  термооптические производные Au (deps/dT_e, deps/dT_l), которые Palpant 2008
  не дает как отдельные константы. Тепловой (thermal-lens) вклад для нс/CW
  берется готовым в виде эквивалентных коэффициентов gamma_th, beta_th
  (Palpant 2008, рис. 2): gamma_th до ~-1.8e-9 см^2/Вт, beta_th до
  ~-0.55e-3 см/Вт у SPR (~2.25 эВ), 7.5 нс, P_abs0=5e17 Вт/м^3 - того же
  порядка, что электронный Керр Ryasnyanskiy (см. Obsidian-заметку
  "Нелинейность наночастиц Au для EMA и расчетов композита", ветка B).

Граница применимости: TTM корректна при tau_p >~ 1 ps; для более коротких
импульсов - athermal-режим (Boltzmann), TTM завышает отклик (Guillet 2009).

Запуск (из корня репозитория):
    python composite_ema/scripts/ttm_hot_electron.py
Выход: composite_ema/results/ttm_hot_electron_{dynamics,regimes}.{csv,png}
       + печать самопроверок.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from composite_ema import RESULTS_DIR

# --- параметры золота (Guillet 2009, первичные) ---------------------------
GAMMA_E = 66.0        # Дж·м^-3·К^-2, C_e = GAMMA_E * T_e
C_L = 2.44e6          # Дж·м^-3·К^-1
T_DEBYE = 170.0       # К
G_EP = 3.0e16         # Вт·м^-3·К^-1 (electron-phonon coupling)
T0 = 300.0            # К, начальная температура


_FWHM_TO_SIGMA = 1.0 / (2.0 * math.sqrt(2.0 * math.log(2.0)))


def p_abs_gaussian(t: np.ndarray, e_abs: float, tau_fwhm: float, t_center: float) -> np.ndarray:
    """Гауссова объемная плотность поглощенной мощности (Вт/м^3), векторно.

    Нормирована так, что интеграл по времени = e_abs (Дж/м^3).
    tau_fwhm - полная ширина на полувысоте (с).
    """
    sigma = tau_fwhm * _FWHM_TO_SIGMA
    norm = e_abs / (sigma * math.sqrt(2.0 * math.pi))
    return norm * np.exp(-0.5 * ((t - t_center) / sigma) ** 2)


def _derivs(te: float, tl: float, p_abs: float) -> tuple[float, float]:
    ce = GAMMA_E * te
    dte = (-G_EP * (te - tl) + p_abs) / ce
    dtl = (G_EP * (te - tl)) / C_L
    return dte, dtl


@dataclass
class TTMResult:
    t: np.ndarray
    te: np.ndarray
    tl: np.ndarray
    e_abs: float
    tau_fwhm: float
    dt: float


def solve_ttm(
    e_abs: float,
    tau_fwhm: float,
    t_end: float | None = None,
    dt: float | None = None,
) -> TTMResult:
    """Интегрирует TTM методом RK4 с фиксированным шагом.

    Шаг выбирается << характерного времени e-ph tau_ep ~ C_e/G ~ 0.6 ps,
    чтобы явный RK4 оставался устойчивым (система жесткая).
    """
    # Явный RK4 для жесткой системы устойчив при dt << tau_ep = C_e/G ~ 0.66 ps;
    # шаг привязан именно к tau_ep (а не к длине импульса), иначе у длинных
    # импульсов dt был бы больше tau_ep и решение разошлось бы.
    tau_ep = GAMMA_E * T0 / G_EP  # ~6.6e-13 с при T0=300 К
    if dt is None:
        dt = tau_ep / 20.0
    if t_end is None:
        t_end = 4.0 * tau_fwhm + 20.0e-12  # захватывает пик и часть релаксации
    t_center = 2.0 * tau_fwhm

    sigma = tau_fwhm * _FWHM_TO_SIGMA
    norm = e_abs / (sigma * math.sqrt(2.0 * math.pi))
    inv2s2 = 0.5 / (sigma * sigma)

    def p_scalar(tt: float) -> float:
        d = tt - t_center
        return norm * math.exp(-inv2s2 * d * d)

    n = int(math.ceil(t_end / dt)) + 1
    t = np.linspace(0.0, (n - 1) * dt, n)
    te = np.empty(n)
    tl = np.empty(n)
    te[0] = T0
    tl[0] = T0

    for i in range(n - 1):
        ti = t[i]
        te_i, tl_i = te[i], tl[i]
        p1 = p_scalar(ti)
        p2 = p_scalar(ti + 0.5 * dt)
        p4 = p_scalar(ti + dt)

        k1e, k1l = _derivs(te_i, tl_i, p1)
        k2e, k2l = _derivs(te_i + 0.5 * dt * k1e, tl_i + 0.5 * dt * k1l, p2)
        k3e, k3l = _derivs(te_i + 0.5 * dt * k2e, tl_i + 0.5 * dt * k2l, p2)
        k4e, k4l = _derivs(te_i + dt * k3e, tl_i + dt * k3l, p4)

        te[i + 1] = te_i + dt / 6.0 * (k1e + 2 * k2e + 2 * k3e + k4e)
        tl[i + 1] = tl_i + dt / 6.0 * (k1l + 2 * k2l + 2 * k3l + k4l)

    return TTMResult(t=t, te=te, tl=tl, e_abs=e_abs, tau_fwhm=tau_fwhm, dt=dt)


def internal_energy_increase(te: float, tl: float) -> float:
    """U(T) - U(T0), Дж/м^3: электроны U_e=(1/2)gamma_e T^2, решетка U_l=C_l T."""
    u_e = 0.5 * GAMMA_E * (te**2 - T0**2)
    u_l = C_L * (tl - T0)
    return u_e + u_l


def eph_relaxation_time(res: TTMResult) -> float:
    """1/e-время спада (T_e - T_l) после окончания импульса."""
    t_center = 2.0 * res.tau_fwhm
    start = t_center + 2.0 * res.tau_fwhm  # заведомо после импульса
    mask = res.t >= start
    diff = res.te[mask] - res.tl[mask]
    t_sub = res.t[mask]
    if len(diff) < 3 or diff[0] <= 0:
        return float("nan")
    target = diff[0] / math.e
    for j in range(1, len(diff)):
        if diff[j] <= target:
            # линейная интерполяция момента пересечения
            frac = (diff[j - 1] - target) / (diff[j - 1] - diff[j])
            return (t_sub[j - 1] + frac * (t_sub[j] - t_sub[j - 1])) - t_sub[0]
    return float("nan")


def self_test() -> None:
    tau_ep_analytic = GAMMA_E * T0 / G_EP
    print(f"[self_test] tau_ep ~ C_e/G (T0={T0:.0f} К) = {tau_ep_analytic*1e12:.3f} ps")

    # Представительный пс-импульс: E_abs подобрана под пик T_e ~ 1000 К.
    res = solve_ttm(e_abs=3.0e7, tau_fwhm=2.0e-12)

    # 1) Сохранение энергии: прирост внутренней энергии баней = поглощено.
    u_final = internal_energy_increase(res.te[-1], res.tl[-1])
    rel_err = abs(u_final - res.e_abs) / res.e_abs
    assert rel_err < 0.01, (u_final, res.e_abs, rel_err)
    print(f"[self_test] OK: сохранение энергии U_fin={u_final:.3e} vs E_abs={res.e_abs:.3e} "
          f"(отклонение {rel_err:.2%}).")

    # 2) Электроны всегда горячее решетки (при нагреве через электроны).
    assert np.all(res.te >= res.tl - 1e-6), "T_e < T_l - нефизично для e-нагрева"
    print(f"[self_test] OK: T_e >= T_l на всей траектории; пик T_e={res.te.max():.0f} К, "
          f"пик T_l={res.tl.max():.0f} К.")

    # 3) e-ph время в пс-диапазоне, согласуется с Hua 2015 (1.8-5.7 ps).
    tau_ep = eph_relaxation_time(res)
    assert 0.2e-12 < tau_ep < 8.0e-12, tau_ep
    print(f"[self_test] OK: 1/e-время (T_e-T_l) = {tau_ep*1e12:.2f} ps "
          f"(Hua 2015 e-ph: 1.8-5.7 ps; аналитика C_e/G ~ 0.7-2 ps).")

    # 4) Квазиравновесие для длинного импульса. Онсет квазиравновесия задается
    #    масштабом C_l/G ~ 81 ps (а НЕ tau_ep=C_e/G~0.66 ps): установившийся
    #    разрыв T_e-T_l ~ P/G, а накопленный перегрев решетки ~ P*tau_p/C_l,
    #    их отношение = C_l/(G*tau_p). При tau_p=2 нс (режим Palpant 7.5 нс)
    #    это ~0.01 -> lattice heating / thermal-lens, а не hot-electron.
    cl_over_g = C_L / G_EP
    print(f"[self_test] онсет квазиравновесия C_l/G = {cl_over_g*1e12:.0f} ps "
          f"(в отличие от tau_ep=C_e/G={tau_ep_analytic*1e12:.2f} ps).")
    res_long = solve_ttm(e_abs=3.0e7, tau_fwhm=2.0e-9)
    heating = res_long.tl.max() - T0
    gap = float(np.max(res_long.te - res_long.tl))
    ratio = gap / heating
    assert ratio < 0.05, (gap, heating, ratio)
    print(f"[self_test] OK: длинный импульс (2 ns): max(T_e-T_l)={gap:.2f} К << "
          f"перегрев решетки {heating:.0f} К (отношение {ratio:.1%}) - "
          f"квазиравновесие, режим thermal-lens (Palpant), не hot-electron.")


def run_dynamics() -> None:
    """Динамика T_e(t), T_l(t) для пс-импульса (hot-electron режим)."""
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    res = solve_ttm(e_abs=3.0e7, tau_fwhm=2.0e-12)

    # прореживаем для CSV/графика (шаг интегрирования очень мелкий)
    stride = max(1, len(res.t) // 2000)
    t_ps = res.t[::stride] * 1e12

    csv_path = RESULTS_DIR / "ttm_hot_electron_dynamics.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["t_ps", "T_e_K", "T_l_K"])
        for tp, te, tl in zip(t_ps, res.te[::stride], res.tl[::stride]):
            writer.writerow([f"{tp:.4f}", f"{te:.2f}", f"{tl:.2f}"])
    print(f"[dynamics] CSV записан: {csv_path}")

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    ax.plot(t_ps, res.te[::stride], color="#c0392b", label=r"$T_e$ (электроны)")
    ax.plot(t_ps, res.tl[::stride], color="#2c3e50", label=r"$T_l$ (решетка)")
    ax.axhline(T0, color="gray", linewidth=0.6, linestyle=":")
    tau_ep = eph_relaxation_time(res)
    ax.set_xlabel("время, пс")
    ax.set_ylabel("температура, К")
    ax.set_title(
        f"TTM Au: импульс 2 пс, E_abs=3e7 Дж/м³ (G={G_EP:.0e}, γ_e={GAMMA_E:.0f}, "
        f"C_l={C_L:.2e}); τ_e-ph≈{tau_ep*1e12:.1f} пс",
        fontsize=9,
    )
    ax.legend()
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    png_path = RESULTS_DIR / "ttm_hot_electron_dynamics.png"
    fig.savefig(png_path, dpi=150)
    print(f"[dynamics] PNG записан: {png_path}")


def run_regime_map() -> None:
    """Скан пиковых T_e, T_l и разрыва по длительности импульса.

    Показывает границу hot-electron (короткие импульсы: T_e >> T_l) ->
    thermal-lens (длинные: T_e ≈ T_l). E_abs фиксирована.
    """
    e_abs = 3.0e7
    taus = np.array(
        [0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0, 100.0, 300.0, 1000.0, 2000.0]
    ) * 1e-12

    rows = []
    for tau in taus:
        res = solve_ttm(e_abs=e_abs, tau_fwhm=float(tau))
        te_peak = float(res.te.max())
        tl_peak = float(res.tl.max())
        gap = float(np.max(res.te - res.tl))
        regime = "hot-electron" if gap > 0.3 * (tl_peak - T0) else "quasi-eq (thermal-lens)"
        ttm_valid = "TTM ok" if tau >= 1.0e-12 else "athermal (TTM завышает)"
        rows.append(
            {
                "tau_fwhm_ps": tau * 1e12,
                "Te_peak_K": te_peak,
                "Tl_peak_K": tl_peak,
                "max_Te_minus_Tl_K": gap,
                "regime": regime,
                "validity": ttm_valid,
            }
        )

    csv_path = RESULTS_DIR / "ttm_hot_electron_regimes.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for r in rows:
            r_fmt = {k: (f"{v:.3f}" if isinstance(v, float) else v) for k, v in r.items()}
            writer.writerow(r_fmt)
    print(f"[regimes] CSV записан: {csv_path}")
    for r in rows:
        print(
            f"[regimes]   tau={r['tau_fwhm_ps']:6.2f} ps: T_e^peak={r['Te_peak_K']:6.0f} К, "
            f"T_l^peak={r['Tl_peak_K']:6.0f} К, max(T_e-T_l)={r['max_Te_minus_Tl_K']:6.0f} К "
            f"-> {r['regime']} [{r['validity']}]"
        )

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    tau_ps = [r["tau_fwhm_ps"] for r in rows]
    ax.semilogx(tau_ps, [r["Te_peak_K"] for r in rows], "o-", color="#c0392b", label=r"$T_e^{peak}$")
    ax.semilogx(tau_ps, [r["Tl_peak_K"] for r in rows], "s-", color="#2c3e50", label=r"$T_l^{peak}$")
    ax.axvline(1.0, color="gray", linestyle="--", linewidth=1)
    ax.text(1.05, ax.get_ylim()[1] * 0.6, "TTM ok →\n← athermal", fontsize=7, color="gray")
    ax.set_xlabel("длительность импульса FWHM, пс")
    ax.set_ylabel("пиковая температура, К")
    ax.set_title(
        "TTM Au: hot-electron (короткие импульсы) → thermal-lens квазиравновесие "
        "(длинные); E_abs=3e7 Дж/м³",
        fontsize=8.5,
    )
    ax.legend()
    ax.grid(True, which="both", alpha=0.25)
    fig.tight_layout()
    png_path = RESULTS_DIR / "ttm_hot_electron_regimes.png"
    fig.savefig(png_path, dpi=150)
    print(f"[regimes] PNG записан: {png_path}")


if __name__ == "__main__":
    self_test()
    run_dynamics()
    run_regime_map()
