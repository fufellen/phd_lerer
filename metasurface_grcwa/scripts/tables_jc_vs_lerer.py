"""Таблицы оптических постоянных металлов программы Лерера против Джонсона-Кристи: насколько от выбора
таблицы зависят результаты статьи о влиянии модели эффективной среды (скелет по докладу SFM-2026).

Основной остаётся таблица Лерера (материальная библиотека VIBR_T: Au = Au.txt, Ag = Ag_.c, Cu = EPS_MET.c,
ниже 517 нм у меди линейная экстраполяция); Джонсон-Кристи считается параллельно, рядом.

Что считается:
- material: на уровне проницаемости композита (матрица n = 1,77) для Au, Ag, Cu и обеих таблиц -
  условие Фрёлиха, максимумы Im eps_eff по правилам смешивания при C = 0,1, полоса Бруггемана при C = 0,1
  и 0,2 (интервал длин волн, где смесь поглощает даже при вещественной eps металла), затухание плоской
  волны в композите и интервалы непрозрачности слоёв 2,1 и 0,2 мкм (10 дБ за один проход);
- bare: столбиковый поглотитель со столбиками и без них (тот же стек без слоя столбиков), Максвелл-Гарнетт,
  обе таблицы золота; без столбиков - grcwa и, для контроля, матрица переноса;
- summary: сводка обоих расчётов плюс готовых CSV (четыре правила на столбиках, тонкая структура статьи,
  материалы) в одну таблицу «Лерер / Джонсон-Кристи»;
- plot: рисунок сравнения.

Запуск: python tables_jc_vs_lerer.py [material|bare|summary|plot|all]
Выход: results/tables_jc_vs_lerer_*.txt|csv|png
"""
from __future__ import annotations

import cmath
import csv
import io
import math
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import lerer_grcwa as L  # noqa: E402

sys.path.insert(0, str(HERE.parents[1] / "vibr_t" / "scripts"))
import vibr_materials as VM  # noqa: E402

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

EPS_H = complex(L.N_HOST ** 2, 0.0)
RES = L.RESULTS
TABLES = ("lerer", "jc")
TABLE_RU = {"lerer": "таблица Лерера", "jc": "Джонсон-Кристи"}

# оптическая конвенция exp(-i w t), Im eps >= 0
METAL = {
    ("Au", "lerer"): L.eps_au_lerer,
    ("Ag", "lerer"): lambda lam: VM.libmat(2, lam).conjugate(),
    ("Cu", "lerer"): lambda lam: VM.libmat(1, lam).conjugate(),
    ("Au", "jc"): L.eps_au_jc,
    ("Ag", "jc"): L.eps_ag_jc,
    ("Cu", "jc"): L.eps_cu_jc,
}
RANGE = {"lerer": (200.0, 2000.0), "jc": (187.9, 1937.0)}
CU_LERER_LO = 517.0

RULES = {
    "MG": lambda ep, f, lam: L.mg(EPS_H, ep, f),
    "Bruggeman": lambda ep, f, lam: L.bruggeman(EPS_H, ep, f),
    "Lerer": lambda ep, f, lam: L.lerer_formula(EPS_H, ep, f),
    "MLWA10": lambda ep, f, lam: L.mlwa_mg(EPS_H, ep, f, lam, 10.0),
    "MLWA15": lambda ep, f, lam: L.mlwa_mg(EPS_H, ep, f, lam, 15.0),
    "MLWA30": lambda ep, f, lam: L.mlwa_mg(EPS_H, ep, f, lam, 30.0),
    "MLWA50": lambda ep, f, lam: L.mlwa_mg(EPS_H, ep, f, lam, 50.0),
}
OPAQUE = {"столбики, 2,1 мкм": 2.1, "решётка из письма, 1,05 мкм": 1.05, "тонкая структура, 0,2 мкм": 0.2}


def atten_db_per_um(eps: complex, lam_nm: float) -> float:
    """Затухание мощности плоской волны в однородной среде, дБ/мкм."""
    return 10.0 * math.log10(math.e) * 4.0 * math.pi * L.n_of(eps).imag / (lam_nm / 1000.0)


def intervals(lam: np.ndarray, mask: np.ndarray) -> str:
    """Интервалы длин волн, где mask истинна, строкой «a-b, c-d»."""
    out, start = [], None
    for x, m in zip(lam, mask):
        if m and start is None:
            start = x
        if not m and start is not None:
            out.append((start, prev))
            start = None
        prev = x
    if start is not None:
        out.append((start, prev))
    return ", ".join("%.0f-%.0f" % ab for ab in out) if out else "нет"


def bruggeman_band_t(C: float) -> tuple[float, float]:
    """Корни (3C-1)^2 t^2 + [2(3C-1)(2-3C)+8] t + (2-3C)^2 = 0, t = eps_p / eps_m: между ними смесь
    Бруггемана поглощает и при вещественной eps_p."""
    a = (3 * C - 1) ** 2
    b = 2 * (3 * C - 1) * (2 - 3 * C) + 8
    c = (2 - 3 * C) ** 2
    d = math.sqrt(b * b - 4 * a * c)
    r = sorted(((-b - d) / (2 * a), (-b + d) / (2 * a)))
    return r[0], r[1]


def run_material():
    lines = ["Проницаемость композита: матрица n = %.2f, наночастицы - Au, Ag, Cu; оптическая конвенция, Im > 0."
             % L.N_HOST,
             "Таблица Лерера = библиотека VIBR_T (Au = Au.txt, Ag = Ag_.c, Cu = EPS_MET.c, ниже 517 нм - "
             "экстраполяция); Джонсон-Кристи = composite_ema.", ""]
    rows = []
    lam = np.arange(400.0, 900.01, 1.0)
    for metal in ("Au", "Ag", "Cu"):
        lines.append("== %s" % metal)
        for tab in TABLES:
            fn = METAL[(metal, tab)]
            ep = np.array([fn(x) for x in lam])
            # условие Фрёлиха: Re eps_p = -2 eps_m
            g = ep.real + 2.0 * EPS_H.real
            j = int(np.where(np.diff(np.sign(g)) != 0)[0][0])
            lam_f = lam[j] - g[j] * (lam[j + 1] - lam[j]) / (g[j + 1] - g[j])
            extra = ""
            if metal == "Cu" and tab == "lerer" and lam_f < CU_LERER_LO:
                extra = " (в области экстраполяции таблицы)"
            lines.append("  %-15s условие Фрёлиха %.1f нм%s" % (TABLE_RU[tab], lam_f, extra))
            for rule in RULES:
                ee = np.array([RULES[rule](e, 0.10, x) for e, x in zip(ep, lam)])
                i = int(np.argmax(ee.imag))
                at = np.array([atten_db_per_um(e, x) for e, x in zip(ee, lam)])
                k = int(np.argmax(at))
                rows.append(dict(metal=metal, table=tab, rule=rule, lam_max_im=lam[i], max_im=ee.imag[i],
                                 lam_max_att=lam[k], max_att=at[k]))
                op = "; ".join("%s: %s" % (name, intervals(lam, at > 10.0 / d)) for name, d in OPAQUE.items())
                lines.append("  %-15s %-10s max Im eps_eff %.3f при %3.0f нм; max затухание %.1f дБ/мкм при %3.0f нм;"
                             " непрозрачен (нм) - %s" % (TABLE_RU[tab], rule, ee.imag[i], lam[i], at[k], lam[k], op))
        lines.append("")
    # полоса Бруггемана по вещественной части eps металла
    lines.append("== Полоса Бруггемана: интервал, где смесь поглощает при Im eps_p = 0 (по Re eps_p таблицы)")
    for C in (0.10, 0.20):
        t_lo, t_hi = bruggeman_band_t(C)
        lines.append("  C = %.2f: t = eps_p/eps_m от %.2f до %.3f; полюс Максвелла-Гарнетта t = %.3f"
                     % (C, t_lo, t_hi, -(2 + C) / (1 - C)))
        for metal in ("Au", "Ag", "Cu"):
            for tab in TABLES:
                lo, hi = RANGE[tab]
                if metal == "Cu" and tab == "lerer":
                    lo = CU_LERER_LO
                lg = np.arange(max(300.0, lo), hi + 0.01, 1.0)
                t = np.array([METAL[(metal, tab)](x).real for x in lg]) / EPS_H.real
                band = (t >= t_lo) & (t <= t_hi)
                s = intervals(lg, band)
                edge = ""
                if band[-1]:
                    edge = " (до конца таблицы, %.0f нм)" % lg[-1]
                if band[0]:
                    edge += " (с начала таблицы, %.0f нм)" % lg[0]
                lines.append("    %s, %-15s %s%s" % (metal, TABLE_RU[tab], s, edge))
    txt = "\n".join(lines) + "\n"
    print(txt)
    io.open(RES / "tables_jc_vs_lerer_material.txt", "w", encoding="utf-8").write(txt)
    with io.open(RES / "tables_jc_vs_lerer_material.csv", "w", encoding="utf-8", newline="\n") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        for r in rows:
            w.writerow({k: ("%.4f" % v if isinstance(v, float) else v) for k, v in r.items()})


BARE_LAMS = [400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0, 750.0]


def run_bare():
    """Столбиковый поглотитель со столбиками и без, Максвелл-Гарнетт, обе таблицы золота."""
    out = io.open(RES / "tables_jc_vs_lerer_bare.csv", "w", encoding="utf-8", newline="\n")
    out.write("table,variant,lambda_nm,R,T,A,R_dB\n")
    lines = ["Столбиковый поглотитель (период 400 нм, столбики 250 x 250 x 100 нм, композит 2000 нм, прослойка "
             "325 нм, подложка 1,45), Максвелл-Гарнетт, C = 0,1; grcwa nG = 221 со столбиками, без столбиков - "
             "grcwa и матрица переноса (ТММ). «Без столбиков» = столбики убраны, композита 2000 нм; в последней "
             "колонке - слой столбиков заменён сплошным композитом, 2100 нм (так считалась прежняя таблица "
             "заметки о модели COMSOL)"]
    for tab in TABLES:
        au = L.AuModel(tab)
        lines.append("== %s" % TABLE_RU[tab])
        lines.append("  длина волны, нм   A со столбиками   A без столбиков (grcwa / ТММ)   R со, дБ   R без, дБ   "
                     "R без / R со   A сплошного слоя 2100 нм (ТММ)")
        for lam in BARE_LAMS:
            layers, eps_sub = L.stack_pillar_absorber(lam, au, mixing=L.mg, **L.ORIGINAL)
            s1 = L.solve(lam, 0.4, layers, eps_sub=eps_sub, nG=221)
            d = dict(L.ORIGINAL)
            d["h_pil_nm"] = 0.0
            layers0, eps_sub0 = L.stack_pillar_absorber(lam, au, mixing=L.mg, **d)
            s0 = L.solve(lam, 0.4, layers0, eps_sub=eps_sub0, nG=21)
            ec = L.mg(EPS_H, au(lam), 0.10)
            Rt, Tt, At = L.tmm_rta(lam, [(L.n_of(ec), 2000.0), (L.N_SUB, 325.0)], 1.0, L.N_SUB)
            A2100 = L.tmm_rta(lam, [(L.n_of(ec), 2100.0)], 1.0, L.N_SUB)[2]
            for var, s in (("pillars", s1), ("bare", s0)):
                out.write("%s,%s,%.0f,%.6f,%.6f,%.6f,%.3f\n" % (tab, var, lam, s.R, s.T, s.A, 10 * math.log10(max(s.R, 1e-12))))
            lines.append("  %5.0f            %.3f             %.3f / %.3f                   %6.1f     %6.1f      %5.1f"
                         "          %.3f"
                         % (lam, s1.A, s0.A, At, 10 * math.log10(s1.R), 10 * math.log10(s0.R), s0.R / s1.R, A2100))
            out.flush()
            print(lines[-1])
    out.close()
    txt = "\n".join(lines) + "\n"
    io.open(RES / "tables_jc_vs_lerer_bare.txt", "w", encoding="utf-8").write(txt)
    print(txt)


def _read(path: Path):
    return list(csv.DictReader(io.open(path, encoding="utf-8")))


def _series(rows, lam_key="lambda_nm", **flt):
    sel = [r for r in rows if all(r[k] == v for k, v in flt.items())]
    lam = np.array([float(r[lam_key]) for r in sel])
    o = np.argsort(lam)
    return lam[o], {k: np.array([float(r[k]) for r in sel])[o] for k in ("A", "R", "T") if k in sel[0]}


def summarize():
    lines = ["Сводка: результаты статьи по таблице Лерера (основная) и по Джонсону-Кристи (параллельно)", ""]
    # четыре правила на столбиках
    pl = {"lerer": _read(RES / "mixing_rules_on_pillars.csv"), "jc": _read(RES / "mixing_rules_on_pillars_jc.csv")}
    lines.append("== Четыре правила на столбиковом поглотителе, Au, исходная конструкция; A при 600/650/700/750 нм")
    spreads = {}
    for tab in TABLES:
        st = []
        for f in ("MG", "Bruggeman", "MLWA15", "Lerer"):
            lam, s = _series(pl[tab], design="original", formula=f)
            a = s["A"]
            st.append(a)
            lines.append("  %-15s %-10s %s; среднее 400-750 %.3f"
                         % (TABLE_RU[tab], f, " ".join("%.3f" % a[list(lam).index(x)] for x in (600, 650, 700, 750)), a.mean()))
        st = np.array(st)
        sp = st.max(0) - st.min(0)
        on = next((lam[j] for j in range(len(lam)) if sp[j] > 0.02), None)
        spreads[tab] = (lam, sp)
        lines.append("  %-15s разброс четырёх: до 650 нм не больше %.3f; больше 0,02 с %.0f нм; при 750 нм %.3f; "
                     "минимум A по всем правилам до 650 нм %.3f"
                     % (TABLE_RU[tab], sp[lam <= 650].max(), on, sp[-1], st[:, lam <= 650].min()))
    lines.append("")
    # тонкая структура
    ar = {"jc": _read(RES / "mixing_rules_on_article.csv"), "lerer": _read(RES / "mixing_rules_on_article_lerer.csv")}
    lines.append("== Тонкая структура статьи (цилиндры 400 x 100 нм, период 500 нм, слой 100 нм), Au")
    for tab in TABLES:
        st = []
        for f in ("MG", "Bruggeman", "MLWA15", "Lerer"):
            lam, s = _series(ar[tab], formula=f)
            a = s["A"]
            st.append(a)
            i = int(np.argmax(a))
            lines.append("  %-15s %-10s max A %.3f при %.0f нм; A(750) %.3f" % (TABLE_RU[tab], f, a[i], lam[i], a[list(lam).index(750)]))
        st = np.array(st)
        sp = st.max(0) - st.min(0)
        i = int(np.argmax(sp))
        lines.append("  %-15s наибольший разброс четырёх %.3f при %.0f нм" % (TABLE_RU[tab], sp[i], lam[i]))
    lines.append("")
    # материалы
    mt = _read(RES / "materials_on_pillars.csv")
    lines.append("== Материал наночастиц на столбиках, Максвелл-Гарнетт: A(700), A(750), наибольшее R во всей полосе")
    for metal in ("Cu", "Ag", "Au"):
        for tab in TABLES:
            lam, s = _series(mt, metal=metal, table=tab, nG="221")
            a, r = s["A"], s["R"]
            k = int(np.argmax(r))
            lines.append("  %s, %-15s A(700) %.3f, A(750) %.3f, min A %.3f при %.0f нм, max R %.1f дБ при %.0f нм"
                         % (metal, TABLE_RU[tab], a[list(lam).index(700)], a[list(lam).index(750)], a.min(),
                            lam[int(np.argmin(a))], 10 * math.log10(r[k]), lam[k]))
    txt = "\n".join(lines) + "\n"
    print(txt)
    io.open(RES / "tables_jc_vs_lerer_summary.txt", "w", encoding="utf-8").write(txt)


def plot():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(2, 2, figsize=(11, 8.2))
    style = {"lerer": "-", "jc": "--"}
    col = {"MG": "k", "Bruggeman": "tab:red", "MLWA15": "tab:green", "Lerer": "tab:blue",
           "Cu": "tab:orange", "Ag": "0.5", "Au": "goldenrod"}
    name = {"MG": "Максвелл-Гарнетт", "Bruggeman": "Бруггеман", "MLWA15": "MLWA, R = 15 нм", "Lerer": "линейная запись"}

    # (а) оптические постоянные золота
    lam = np.arange(400.0, 900.01, 2.0)
    axn = ax[0, 0].twinx()
    for tab in TABLES:
        nk = np.array([L.n_of(METAL[("Au", tab)](x)) for x in lam])
        ax[0, 0].plot(lam, nk.imag, style[tab], color="goldenrod", label="k, " + TABLE_RU[tab])
        axn.plot(lam, nk.real, style[tab], color="0.35", label="n, " + TABLE_RU[tab])
    ax[0, 0].axvline(650, color="0.7", lw=0.8, ls=":")
    ax[0, 0].set_title("(а) золото: показатели преломления n и поглощения k")
    ax[0, 0].set_xlabel("Длина волны, нм")
    ax[0, 0].set_ylabel("Показатель поглощения k")
    axn.set_ylabel("Показатель преломления n")
    h1, l1 = ax[0, 0].get_legend_handles_labels()
    h2, l2 = axn.get_legend_handles_labels()
    ax[0, 0].legend(h1 + h2, l1 + l2, fontsize=8, loc="upper center")

    # (б) четыре правила на столбиках
    pl = {"lerer": _read(RES / "mixing_rules_on_pillars.csv"), "jc": _read(RES / "mixing_rules_on_pillars_jc.csv")}
    for tab in TABLES:
        for f in ("MG", "Bruggeman", "MLWA15", "Lerer"):
            x, s = _series(pl[tab], design="original", formula=f)
            ax[0, 1].plot(x, s["A"], style[tab], color=col[f], label=name[f] if tab == "lerer" else None)
    ax[0, 1].plot([], [], "-", color="0.3", label=TABLE_RU["lerer"])
    ax[0, 1].plot([], [], "--", color="0.3", label=TABLE_RU["jc"])
    ax[0, 1].set_title("(б) столбиковый поглотитель, Au, четыре правила")
    ax[0, 1].set_xlabel("Длина волны, нм")
    ax[0, 1].set_ylabel("Поглощение")
    ax[0, 1].legend(fontsize=8, loc="lower left")

    # (в) тонкая структура
    ar = {"jc": _read(RES / "mixing_rules_on_article.csv"), "lerer": _read(RES / "mixing_rules_on_article_lerer.csv")}
    for tab in TABLES:
        for f in ("MG", "Bruggeman", "MLWA15", "Lerer"):
            x, s = _series(ar[tab], formula=f)
            m = x <= 900
            ax[1, 0].plot(x[m], s["A"][m], style[tab], color=col[f], label=name[f] if tab == "lerer" else None)
    ax[1, 0].plot([], [], "-", color="0.3", label=TABLE_RU["lerer"])
    ax[1, 0].plot([], [], "--", color="0.3", label=TABLE_RU["jc"])
    ax[1, 0].set_title("(в) тонкая структура статьи, Au, четыре правила")
    ax[1, 0].set_xlabel("Длина волны, нм")
    ax[1, 0].set_ylabel("Поглощение")
    ax[1, 0].set_ylim(0, 1.32)
    ax[1, 0].set_yticks(np.arange(0, 1.01, 0.2))
    ax[1, 0].legend(fontsize=8, loc="upper center", ncol=3)

    # (г) материалы
    mt = _read(RES / "materials_on_pillars.csv")
    for metal in ("Cu", "Ag", "Au"):
        for tab in TABLES:
            x, s = _series(mt, metal=metal, table=tab, nG="221")
            ax[1, 1].plot(x, s["A"], style[tab], color=col[metal], label=metal if tab == "lerer" else None)
    ax[1, 1].plot([], [], "-", color="0.3", label=TABLE_RU["lerer"])
    ax[1, 1].plot([], [], "--", color="0.3", label=TABLE_RU["jc"])
    ax[1, 1].set_title("(г) столбики, Максвелл-Гарнетт, Cu, Ag, Au")
    ax[1, 1].set_xlabel("Длина волны, нм")
    ax[1, 1].set_ylabel("Поглощение")
    ax[1, 1].legend(fontsize=8, loc="lower left")
    for a in ax.flat:
        a.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(RES / "tables_jc_vs_lerer.png", dpi=150)
    print("рисунок:", RES / "tables_jc_vs_lerer.png")


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    t0 = time.time()
    if which in ("material", "all"):
        run_material()
    if which in ("bare", "all"):
        run_bare()
    if which in ("summary", "all"):
        summarize()
    if which in ("plot", "all"):
        plot()
    print("готово за %.0f с" % (time.time() - t0))
