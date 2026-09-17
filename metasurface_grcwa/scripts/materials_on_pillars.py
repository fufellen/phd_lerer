"""Материал наночастиц на столбиковом поглотителе Лерера: P и R для Cu, Ag, Au по формуле «КМ».

Проверка рисунков «Сравнение разных материалов наночастиц. Потери / Коэффициент отражения» из задания
А.М. Лерера к докладу SFM-2026 (01.09.2026). Конструкция исходная (lerer_grcwa.ORIGINAL): период 400 нм,
столбики 250 x 250 x 100 нм, слой композита 2000 нм, прослойка 325 нм, подложка n = 1,45; композит -
наночастицы 10 % в матрице n = 1,77, правило смешивания - полный Максвелл-Гарнетт (у Лерера - «формула
Клаузиуса-Моссотти», КМ; для этой конструкции подпись «КМ» = полный МГ установлена в
mixing_rules_on_pillars.py: 0,566 при 750 нм против 0,57 на рисунке).

Металлы - две базы:
- «lerer»: материальная библиотека программы VIBR_T (порт vibr_t/scripts/vibr_materials.py): Cu = libmat(1),
  таблица EPS_MET.c 517-3900 нм, ниже 517 нм линейная экстраполяция, как в C-коде; Ag = libmat(2), плотная
  таблица Ag_.c; Au = libmat(3), плотная таблица Au_.c (совпадает с Au.txt, которой сверен рисунок четырёх формул);
- «jc»: Johnson & Christy 1972 из composite_ema - для оценки того, насколько вывод зависит от таблицы.

Запуск:
  python materials_on_pillars.py            -> results/materials_on_pillars.csv (дозапись, продолжение прерванного)
  python materials_on_pillars.py --conv     -> сходимость по числу гармоник в контрольных точках
  python materials_on_pillars.py --summary  -> results/materials_on_pillars_summary.txt по готовому CSV
"""
from __future__ import annotations

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

CSV = L.RESULTS / "materials_on_pillars.csv"
SUMMARY = L.RESULTS / "materials_on_pillars_summary.txt"
HEADER = "metal,table,nG,lambda_nm,Re_eps_metal,Im_eps_metal,Re_eps_eff,Im_eps_eff,R,T,A,R_dB\n"
EPS_H = complex(L.N_HOST ** 2, 0.0)
NG = 221

# оптическая конвенция exp(-i w t), Im eps >= 0
METALS = {
    ("Cu", "lerer"): lambda lam: VM.libmat(1, lam).conjugate(),
    ("Ag", "lerer"): lambda lam: VM.libmat(2, lam).conjugate(),
    ("Au", "lerer"): lambda lam: VM.libmat(3, lam).conjugate(),
    ("Cu", "jc"): L.eps_cu_jc,
    ("Ag", "jc"): L.eps_ag_jc,
    ("Au", "jc"): L.eps_au_jc,
}

COARSE = [float(x) for x in range(400, 751, 5)]
# уточнение около особенностей, прочитанных на рисунках Лерера
FINE = {
    ("Ag", "lerer"): [float(x) for x in range(446, 471, 1)] + [float(x) for x in range(520, 536, 1)],
    ("Cu", "lerer"): [float(x) for x in range(406, 436, 1)] + [float(x) for x in range(722, 742, 1)],
    ("Au", "lerer"): [float(x) for x in range(530, 551, 1)] + [float(x) for x in range(728, 746, 1)],
}


def solve_point(metal_fn, lam: float, nG: int):
    em = metal_fn(lam)
    layers, eps_sub = L.stack_pillar_absorber(lam, metal_fn, mixing=L.mg, **L.ORIGINAL)
    s = L.solve(lam, L.ORIGINAL["period_nm"] / 1000.0, layers, eps_sub=eps_sub, nG=nG)
    ec = L.mg(EPS_H, em, 0.10)
    return em, ec, s


def done_keys():
    keys = set()
    if CSV.exists():
        for line in io.open(CSV, encoding="utf-8"):
            p = line.strip().split(",")
            if len(p) >= 4 and p[0] != "metal":
                keys.add((p[0], p[1], int(p[2]), float(p[3])))
    return keys


def run(points, nG: int):
    new = not CSV.exists()
    have = done_keys()
    t0 = time.time()
    n_new = 0
    with io.open(CSV, "a", encoding="utf-8", newline="\n") as out:
        if new:
            out.write(HEADER)
        for (metal, table), lams in points:
            fn = METALS[(metal, table)]
            for lam in lams:
                if (metal, table, nG, lam) in have:
                    continue
                em, ec, s = solve_point(fn, lam, nG)
                r_db = 10.0 * math.log10(max(s.R, 1e-12))
                out.write("%s,%s,%d,%.1f,%.6f,%.6f,%.6f,%.6f,%.7f,%.7f,%.7f,%.3f\n"
                          % (metal, table, nG, lam, em.real, em.imag, ec.real, ec.imag, s.R, s.T, s.A, r_db))
                out.flush()
                have.add((metal, table, nG, lam))
                n_new += 1
            print("%s/%s nG=%d готово, всего новых точек %d, %.0f с" % (metal, table, nG, n_new, time.time() - t0),
                  flush=True)


def load(nG: int = NG):
    data = {}
    for line in io.open(CSV, encoding="utf-8"):
        p = line.strip().split(",")
        if p[0] == "metal" or int(p[2]) != nG:
            continue
        key = (p[0], p[1])
        data.setdefault(key, {})[float(p[3])] = (float(p[8]), float(p[9]), float(p[10]), float(p[11]))
    out = {}
    for key, d in data.items():
        lams = np.array(sorted(d))
        arr = np.array([d[l] for l in lams])
        out[key] = (lams, arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3])
    return out


def summary():
    data = load()
    lines = ["Столбиковый поглотитель, формула КМ (полный Максвелл-Гарнетт), C = 10 %%, grcwa nG = %d" % NG, ""]
    for key in sorted(data):
        lams, R, T, A, Rdb = data[key]
        at = lambda arr, l: float(arr[np.argmin(np.abs(lams - l))])
        i_rmax = int(np.argmax(Rdb))
        blue = lams <= 500
        i_rmax_blue = int(np.flatnonzero(blue)[np.argmax(Rdb[blue])])
        lo = (lams >= 430) & (lams <= 500)
        i_amin_blue = int(np.flatnonzero(lo)[np.argmin(A[lo])])
        red = lams >= 550
        lines.append("%s/%s: P(700) = %.3f, P(750) = %.3f, среднее P 400-750 = %.3f" % (key[0], key[1], at(A, 700.0),
                     at(A, 750.0), float(np.mean(A[np.isin(lams, np.arange(400, 751, 5))]))))
        lines.append("   min P на 430-500 нм: %.3f при %.0f нм; max R_dB на 400-500 нм: %.2f дБ при %.0f нм"
                     % (A[i_amin_blue], lams[i_amin_blue], Rdb[i_rmax_blue], lams[i_rmax_blue]))
        lines.append("   max R_dB на всей полосе: %.2f дБ при %.0f нм (R = %.4f); max R_dB на 550-750 нм: %.2f дБ при %.0f нм"
                     % (Rdb[i_rmax], lams[i_rmax], R[i_rmax], float(np.max(Rdb[red])),
                        lams[red][int(np.argmax(Rdb[red]))]))
        above = lams[(lams <= 520) & (Rdb > -28.0)]
        lines.append("   R_dB > -28 дБ на 400-520 нм при: %s" % (", ".join("%.0f" % l for l in above) if above.size else "нигде"))
        if key == ("Au", "lerer"):
            mid = (lams >= 520) & (lams <= 560)
            j = int(np.flatnonzero(mid)[np.argmax(Rdb[mid])])
            lines.append("   локальный max R_dB на 520-560 нм: %.2f дБ при %.0f нм" % (Rdb[j], lams[j]))
    conv = {}
    for line in io.open(CSV, encoding="utf-8"):
        p = line.strip().split(",")
        if p[0] != "metal":
            conv.setdefault((p[0], p[1], float(p[3])), {})[int(p[2])] = (float(p[10]), float(p[11]))
    lines += ["", "Сходимость по числу гармоник (P; R, дБ):"]
    for (metal, table, lam), d in sorted(conv.items()):
        if len(d) > 1:
            lines.append("   %s/%s %.0f нм: %s" % (metal, table, lam, "; ".join(
                "nG=%d: %.4f, %.2f" % (g, d[g][0], d[g][1]) for g in sorted(d))))
    txt = "\n".join(lines) + "\n"
    io.open(SUMMARY, "w", encoding="utf-8", newline="\n").write(txt)
    print(txt)


def convergence():
    pts = [("Cu", "lerer", 420.0), ("Ag", "lerer", 458.0), ("Au", "lerer", 540.0), ("Cu", "lerer", 700.0),
           ("Ag", "lerer", 750.0)]
    for nG in (121, 221, 361):
        for metal, table, lam in pts:
            run([((metal, table), [lam])], nG)


if __name__ == "__main__":
    t_start = time.time()
    if "--summary" in sys.argv:
        summary()
    elif "--conv" in sys.argv:
        convergence()
    elif "--time" in sys.argv:
        for nG in (121, 221):
            t = time.time()
            em, ec, s = solve_point(METALS[("Ag", "lerer")], 458.0, nG)
            print("nG=%d: R=%.5f T=%.5f A=%.5f, %.1f с" % (nG, s.R, s.T, s.A, time.time() - t), flush=True)
        a = np.array([abs(VM.libmat(3, l).conjugate() - L.eps_au_lerer(l)) for l in range(400, 751)])
        print("Au_.c против Au.txt, 400-750 нм: max |d eps| = %.2e" % a.max())
    else:
        order = [("Ag", "lerer"), ("Cu", "lerer"), ("Au", "lerer"), ("Ag", "jc"), ("Cu", "jc"), ("Au", "jc")]
        run([(k, COARSE) for k in order], NG)
        run([(k, FINE[k]) for k in FINE], NG)
        summary()
    print("всего %.0f с" % (time.time() - t_start))
