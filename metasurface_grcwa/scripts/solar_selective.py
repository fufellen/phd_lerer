"""Второе применение столбикового поглотителя: спектрально-селективный поглотитель солнечного излучения.

«Мои предложения», пункт 3: у Лерера оптические поглотители, «быть может найти и другое применение».
Композит золотых наночастиц в диэлектрике - кермет, а поглотитель с золотым зеркалом снизу - ровно
устройство селективных солнечных покрытий: поглощать солнечный спектр (0,3...2,5 мкм) и не излучать в
тепловом ИК (2,5...50 мкм). Здесь для трёх конструкций считаются:
- спектр поглощения A(lambda) на 0,3...4 мкм (grcwa) и 2,5...50 мкм (grcwa, малое число гармоник);
- солнечная поглощательная способность alpha_s по спектру AM1.5G (ASTM G173, из pvlib);
- тепловая излучательная способность eps_th(T) по закону Кирхгофа (eps = A) с весом Планка при
  373, 473, 573 и 773 К, нормальное падение;
- к.п.д. фототермического преобразования eta = alpha_s - eps_th sigma T^4 / (C q_sun).
Конструкции: исходная (подложка n = 1,45), предложенная (два слоя композита, золотое зеркало) и
предложенная без зеркала (два слоя композита на подложке).

Допущения, которые нужно помнить при чтении чисел: матрица n = 1,77 и прослойка n = 1,45 без потерь во
всём диапазоне (фононное поглощение оксидов в ИК не учтено); золото за пределами таблицы (2 мкм) - модель
Друде, подобранная по последним 500 нм таблицы, чувствительность к ней проверена моделью Ordal 1985;
для конструкций на подложке в ИК излучательная способность взята как 1 - R (толстая стеклянная подложка
непрозрачна в тепловом ИК), а в солнечной полосе поглощение - как 1 - R - T (прозрачная подложка).

Запуск: python solar_selective.py
Выход: results/solar_selective_spectra.csv, results/solar_selective_summary.txt
"""
from __future__ import annotations

import io
import math
import sys
import time

import numpy as np

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import lerer_grcwa as L  # noqa: E402

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

SIGMA = 5.670374419e-8
H_PLANCK, K_B = 6.62607015e-34, 1.380649e-23
Q_SUN = 1000.0

DESIGNS = {
    "original": (L.ORIGINAL, False),
    "improved": (L.IMPROVED, True),
    "improved_no_mirror": (L.IMPROVED_NO_MIRROR, False),
}
DESIGN_RU = {
    "original": "исходная конструкция, подложка",
    "improved": "предложенная конструкция, золотое зеркало",
    "improved_no_mirror": "предложенная без зеркала, подложка",
}
LAM_SOLAR = sorted(set(list(range(300, 800, 10)) + list(range(800, 2000, 20)) + list(range(2000, 4001, 50))))
LAM_IR = np.geomspace(2500.0, 50000.0, 40)
TEMPS = (373.0, 473.0, 573.0, 773.0)


def n_g(period_nm: float, lam_nm: float) -> int:
    ratio = lam_nm / period_nm
    if ratio < 2.5:
        return 221 if period_nm >= 400 else 121
    if ratio < 6:
        return 121 if period_nm >= 400 else 61
    return 25


def spectrum(design: str, au, lams):
    d, mirror = DESIGNS[design]
    rows = []
    for lam in lams:
        layers, eps_sub = L.stack_pillar_absorber(float(lam), au, **d)
        s = L.solve(float(lam), d["period_nm"] / 1000.0, layers, eps_sub=eps_sub, nG=n_g(d["period_nm"], float(lam)))
        rows.append((float(lam), s.R, s.T, s.A, 1.0 - s.R))
    return np.array(rows)


def am15g():
    from pvlib import spectrum as ps
    ref = ps.get_reference_spectra()
    return ref.index.to_numpy(dtype=float), ref["global"].to_numpy(dtype=float)


def planck(lam_nm, T):
    lam = lam_nm * 1e-9
    x = H_PLANCK * 2.99792458e8 / (lam * K_B * T)
    return 1.0 / (lam ** 5 * np.expm1(x))


def weighted(lam_nm, a, w):
    return float(np.trapezoid(a * w, lam_nm) / np.trapezoid(w, lam_nm))


def solar_absorptance(lam_nm, a):
    wl, g = am15g()
    sel = (wl >= lam_nm.min()) & (wl <= lam_nm.max())
    ai = np.interp(wl[sel], lam_nm, a)
    return weighted(wl[sel], ai, g[sel]), float(np.trapezoid(g[sel], wl[sel]) / np.trapezoid(g, wl))


def emittance(lam_nm, a, T):
    grid = np.geomspace(300.0, 50000.0, 3000)
    ai = np.interp(grid, lam_nm, a)
    return weighted(grid, ai, planck(grid, T))


def main():
    t0 = time.time()
    au = L.AuModel("lerer", "fit")
    au_ord = L.AuModel("lerer", "ordal")
    au_jc = L.AuModel("jc", "fit")
    print(au.describe()); print(au_ord.describe()); print(au_jc.describe())
    # проверка числа гармоник в ИК
    d, _ = DESIGNS["improved"]
    layers, eps_sub = L.stack_pillar_absorber(2500.0, au, **d)
    for nG in (25, 61):
        s = L.solve(2500.0, 0.2, layers, eps_sub=eps_sub, nG=nG)
        print("ИК-проверка 2,5 мкм, предложенная конструкция, nG=%d: R=%.5f" % (nG, s.R))
    out = io.open(L.RESULTS / "solar_selective_spectra.csv", "w", encoding="utf-8", newline="\n")
    out.write("design,au_model,lambda_nm,R,T,A_layers,A_total\n")
    lines = []
    for design in DESIGNS:
        mirror = DESIGNS[design][1]
        sol = spectrum(design, au, LAM_SOLAR)
        ir = spectrum(design, au, LAM_IR)
        ir_ord = spectrum(design, au_ord, LAM_IR)
        for tag, arr in (("lerer_fit", sol), ("lerer_fit", ir), ("lerer_ordal", ir_ord)):
            for lam, R, T, A, At in arr:
                out.write("%s,%s,%.1f,%.6f,%.6f,%.6f,%.6f\n" % (design, tag, lam, R, T, A, At))
        out.flush()
        # солнечная полоса: с зеркалом 1 - R, на подложке 1 - R - T
        a_sol = sol[:, 4] if mirror else sol[:, 3]
        alpha, cover = solar_absorptance(sol[:, 0], a_sol)
        # ИК: eps = 1 - R в обоих случаях (толстая подложка непрозрачна в тепловом ИК)
        lam_all = np.concatenate([sol[:, 0], ir[:, 0]])
        a_all = np.concatenate([sol[:, 4], ir[:, 4]])
        a_all_ord = np.concatenate([sol[:, 4], ir_ord[:, 4]])
        # вариант «подложка прозрачна и сзади ничего не поглощает»: A = 1 - R - T
        a_all_layers = np.concatenate([sol[:, 3], ir[:, 3]])
        lines.append("%s (%s)" % (design, DESIGN_RU[design]))
        lines.append("  alpha_s (AM1.5G, %.1f %% спектра в окне %d...%d нм) = %.4f" % (100 * cover, LAM_SOLAR[0], LAM_SOLAR[-1], alpha))
        lines.append("  среднее A на 400...750 = %.4f, на 750...1500 = %.4f, на 1500...2500 = %.4f"
                     % (np.mean(a_sol[(sol[:, 0] >= 400) & (sol[:, 0] <= 750)]),
                        np.mean(a_sol[(sol[:, 0] >= 750) & (sol[:, 0] <= 1500)]),
                        np.mean(a_sol[(sol[:, 0] >= 1500) & (sol[:, 0] <= 2500)])))
        for T in TEMPS:
            e = emittance(lam_all, a_all, T)
            e_ord = emittance(lam_all, a_all_ord, T)
            e_lay = emittance(lam_all, a_all_layers, T)
            eta1 = alpha - e * SIGMA * T ** 4 / Q_SUN
            eta10 = alpha - e * SIGMA * T ** 4 / (10 * Q_SUN)
            lines.append("  T = %3.0f K: eps_th = %.4f (Ordal: %.4f; только слои 1-R-T: %.4f); "
                         "eta(1 солнце) = %+.3f, eta(10 солнц) = %+.3f, alpha/eps = %.1f"
                         % (T, e, e_ord, e_lay, eta1, eta10, alpha / e))
        print("\n".join(lines[-6:]))
        print("  [%.0f с]" % (time.time() - t0))
    # чувствительность alpha_s к таблице золота на солнечной полосе (предложенная конструкция)
    sol_jc = spectrum("improved", au_jc, LAM_SOLAR)
    alpha_jc, _ = solar_absorptance(sol_jc[:, 0], sol_jc[:, 4])
    for lam, R, T, A, At in sol_jc:
        out.write("improved,jc_fit,%.1f,%.6f,%.6f,%.6f,%.6f\n" % (lam, R, T, A, At))
    lines.append("предложенная конструкция по Джонсону-Кристи: alpha_s = %.4f" % alpha_jc)
    out.close()
    txt = "\n".join(lines) + "\n"
    print(txt)
    io.open(L.RESULTS / "solar_selective_summary.txt", "w", encoding="utf-8").write(txt)
    print("готово за %.0f с" % (time.time() - t0))


if __name__ == "__main__":
    main()
