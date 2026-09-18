"""Правила смешивания композита на метаповерхности статьи: поглощение, отражение и прохождение.

Повод - доклад по теме «Влияние модели эффективной среды на расчётные характеристики поглощающих
метаповерхностей с плазмонными наночастицами». В рукописи этой статьи влияние модели эффективной среды
показано только на уровне самой проницаемости (максимумы Im eps_eff), а спектры решётки посчитаны одной
формулой Максвелла-Гарнетта. Здесь та же решётка статьи посчитана при каждом правиле смешивания, то есть
разброс моделей доведён до расчётной характеристики устройства.

Геометрия статьи: подложка n = 1,45, сплошной слой композита 100 нм, над ним квадратная решётка цилиндров
того же композита высотой 100 нм, период 500 нм, диаметр 400 нм, сверху воздух; композит - Au 10 % в
матрице n = 1,77; золото по Джонсону-Кристи (таблица composite_ema), нормальное падение.

Правила: MG - полный Максвелл-Гарнетт (в коде Лерера названный формулой Клаузиуса-Моссотти), Bruggeman -
симметричная модель Бруггемана, MLWA10 и MLWA15 - Максвелл-Гарнетт через поляризуемость сферы с
поправкой MLWA при радиусе 10 и 15 нм, Lerer - линейная по концентрации запись Лерера.

Запуск: python mixing_rules_on_article.py [conv|main|summary|all]
Выход: results/mixing_rules_on_article*.csv|txt
"""
from __future__ import annotations

import csv
import io
import sys
import time

import numpy as np

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import lerer_grcwa as L  # noqa: E402

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

PERIOD, DIAM, H_CYL, H_FILM, SPACER = 500.0, 400.0, 100.0, 100.0, 0.0
EPS_H = complex(L.N_HOST ** 2, 0.0)
NG = 121
LAMS = list(range(400, 901, 5))
FORMULAS = ["MG", "Bruggeman", "MLWA10", "MLWA15", "Lerer"]
MIX = dict(L.MIXING)
MIX["MLWA10"] = lambda eh, ep, f, lam: L.mlwa_mg(eh, ep, f, lam, 10.0)
AU = L.AuModel("jc")
OUT = L.RESULTS / "mixing_rules_on_article.csv"

# Записанные ранее значения той же структуры при формуле MG - контроль, что решается та же задача:
# собственная двумерная RCWA (N = 5, 121 гармоника) и сторонний grcwa главы 4 диссертации (11 гармоник).
RECORDED_OWN_RCWA = {550: 0.9235}
RECORDED_GRCWA_CH4 = {500: 0.7226, 530: 0.8789, 550: 0.9175, 570: 0.9070, 590: 0.8633, 620: 0.5321}


def solve_article(lam, eps_c, nG=NG):
    layers, eps_sub = L.stack_cylinder_metasurface(lam, eps_c, eps_c, PERIOD, DIAM, H_CYL, H_FILM, SPACER)
    s = L.solve(lam, PERIOD / 1000.0, layers, eps_sub=eps_sub, nG=nG)
    return s.R, s.T, s.A


def run_conv():
    lines = ["сходимость по числу гармоник nG; A при 450/550/650/750/850 нм"]
    for fname in ("MG", "Bruggeman"):
        for nG in (61, 121, 221, 361):
            t0 = time.time()
            row = []
            for lam in (450, 550, 650, 750, 850):
                eps_c = MIX[fname](EPS_H, AU(lam), L.C_FILL, lam)
                row.append("%.4f" % solve_article(lam, eps_c, nG)[2])
            lines.append("%-9s nG = %3d (%.0f с): %s" % (fname, nG, time.time() - t0, " ".join(row)))
            print(lines[-1])
    lines.append("контроль MG против записанных значений (nG = %d):" % NG)
    for lam in sorted(set(RECORDED_GRCWA_CH4) | set(RECORDED_OWN_RCWA)):
        eps_c = MIX["MG"](EPS_H, AU(lam), L.C_FILL, lam)
        a = solve_article(lam, eps_c)[2]
        rec = []
        if lam in RECORDED_OWN_RCWA:
            rec.append("своя RCWA %.4f" % RECORDED_OWN_RCWA[lam])
        if lam in RECORDED_GRCWA_CH4:
            rec.append("grcwa главы 4 %.4f" % RECORDED_GRCWA_CH4[lam])
        lines.append("  %d нм: A = %.4f; записано: %s" % (lam, a, ", ".join(rec)))
        print(lines[-1])
    io.open(L.RESULTS / "mixing_rules_on_article_convergence.txt", "w", encoding="utf-8").write("\n".join(lines) + "\n")


def run_main():
    out = io.open(OUT, "w", encoding="utf-8", newline="\n")
    out.write("formula,lambda_nm,Re_eps_eff,Im_eps_eff,R,T,A\n")
    for fname in FORMULAS:
        t0 = time.time()
        for lam in LAMS:
            eps_c = MIX[fname](EPS_H, AU(lam), L.C_FILL, lam)
            R, T, A = solve_article(lam, eps_c)
            out.write("%s,%d,%.6f,%.6f,%.6f,%.6f,%.6f\n" % (fname, lam, eps_c.real, eps_c.imag, R, T, A))
        out.flush()
        print("%-10s готово за %.0f с" % (fname, time.time() - t0))
    out.close()


def summarize():
    rows = list(csv.DictReader(io.open(OUT, encoding="utf-8")))
    lam = np.array(LAMS, dtype=float)
    get = lambda f, key: np.array([float(r[key]) for r in rows if r["formula"] == f])  # noqa: E731
    A = {f: get(f, "A") for f in FORMULAS}
    R = {f: get(f, "R") for f in FORMULAS}
    T = {f: get(f, "T") for f in FORMULAS}
    im = {f: get(f, "Im_eps_eff") for f in FORMULAS}
    lines = ["метаповерхность статьи: период 500 нм, цилиндры 400 x 100 нм, слой 100 нм, подложка 1,45; "
             "Au по Джонсону-Кристи, C = 0,10, n_h = 1,77; nG = %d" % NG]
    at = (450, 500, 550, 600, 650, 700, 750, 800, 850, 900)
    lines.append("формула    макс. A (нм)     среднее 400-900  среднее 400-650  " + " ".join("A(%d)" % x for x in at))
    for f in FORMULAS:
        a = A[f]
        i = int(np.argmax(a))
        lines.append("%-10s %.4f (%3d)    %.4f           %.4f           %s"
                     % (f, a[i], LAMS[i], a.mean(), a[lam <= 650].mean(),
                        " ".join("%.3f " % a[LAMS.index(x)] for x in at)))
    lines.append("максимум Im eps_eff на сетке расчёта (оптическая конвенция, Im > 0):")
    for f in FORMULAS:
        j = int(np.argmax(im[f]))
        lines.append("  %-10s %.3f при %d нм" % (f, im[f][j], LAMS[j]))
    stack = np.array([A[f] for f in FORMULAS])
    main4 = np.array([A[f] for f in ("MG", "Bruggeman", "MLWA15", "Lerer")])
    for name, st in (("все пять", stack), ("четыре (MG, Bruggeman, MLWA15, Lerer)", main4)):
        spread = st.max(axis=0) - st.min(axis=0)
        i = int(np.argmax(spread))
        lines.append("разброс A, %s: наибольший %.3f при %d нм; средний на 400-650 нм %.3f, на 650-900 нм %.3f"
                     % (name, spread[i], LAMS[i], spread[lam <= 650].mean(), spread[lam >= 650].mean()))
        for thr in (0.02, 0.05, 0.10):
            where = [LAMS[j] for j in range(len(LAMS)) if spread[j] > thr]
            lines.append("  разброс больше %.2f: %s" % (thr, ("с %d по %d нм, точек %d" % (where[0], where[-1], len(where)))
                                                            if where else "нигде"))
    for f in FORMULAS[1:]:
        d = A[f] - A["MG"]
        j = int(np.argmax(np.abs(d)))
        lines.append("%-10s - MG: наибольшее отличие A %+.3f при %d нм" % (f, d[j], LAMS[j]))
    lines.append("отражение и прохождение в характерных точках (R / T):")
    for x in (550, 650, 750, 850):
        k = LAMS.index(x)
        lines.append("  %d нм: %s" % (x, "; ".join("%s %.3f / %.3f" % (f, R[f][k], T[f][k]) for f in FORMULAS)))
    txt = "\n".join(lines) + "\n"
    print(txt)
    io.open(L.RESULTS / "mixing_rules_on_article_summary.txt", "w", encoding="utf-8").write(txt)


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    t0 = time.time()
    if which in ("conv", "all"):
        run_conv()
    if which in ("main", "all"):
        run_main()
    if which in ("main", "summary", "all"):
        summarize()
    print("готово за %.0f с" % (time.time() - t0))
