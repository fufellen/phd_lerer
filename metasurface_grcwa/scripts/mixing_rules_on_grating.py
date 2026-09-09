"""Правила смешивания композита на самой метаповерхности: поглощение P(lambda) решётки цилиндров.

Повод - письмо А.М. Лерера от 19.08.2026 «Сравнение 4 формул расчета композита»: его график P(lambda)
для четырёх формул («КМ», «Бруггемана», «Максвелл-Гарнетт R=15 нм», «Лерер») до сих пор был известен
только по растру - число 0,5 у 700...730 нм снято с рисунка. Здесь та же структура посчитана открытым
решателем grcwa (метод фурье-мод).

Геометрия - входной файл того письма: период 500 x 500 нм, цилиндры диаметром 300 нм и высотой 100 нм в
воздухе, под ними слой композита 950 нм, прослойка 325 нм n = 1,45, подложка n = 1,45; композит -
Au 10 % в матрице n = 1,77; нормальное падение, полоса 400...750 нм.

Формулы: MG - полный Максвелл-Гарнетт (в коде Лерера назван формулой Клаузиуса-Моссотти), Bruggeman,
MLWA15 - Максвелл-Гарнетт через поляризуемость сферы с поправкой MLWA при R = 15 нм, Lerer - формула
Лерера; для справки dilute (первый порядок по концентрации) и LL (Лоренц-Лоренц). Две таблицы золота:
Лерера (Au.txt) и Джонсона-Кристи. Отдельно - два сорта наночастиц (исправленная COMPOSITE_3) и
кривая с ошибкой исходной COMPOSITE_3 при одном сорте (расчёт Яковлева).

Запуск: python mixing_rules_on_grating.py [conv|main|extra|all]
Выход: results/mixing_rules_on_grating*.csv|txt
"""
from __future__ import annotations

import io
import sys
import time

import numpy as np

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import lerer_grcwa as L  # noqa: E402

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

sys.path.insert(0, str(L.REPO / "vibr_t" / "scripts"))
try:
    import vibr_composite as vc  # noqa: E402
except Exception:  # pragma: no cover
    vc = None

PERIOD, DIAM, H_CYL, H_FILM, SPACER = 500.0, 300.0, 100.0, 950.0, 325.0
EPS_H = complex(L.N_HOST ** 2, 0.0)
NG = 221
LAMS = list(range(400, 751, 5))
LAMS10 = list(range(400, 751, 10))
FORMULAS = ["MG", "Bruggeman", "MLWA15", "Lerer"]
EXTRA_FORMULAS = ["dilute", "LL"]
AU = {"lerer": L.AuModel("lerer"), "jc": L.AuModel("jc")}


def solve_grating(lam, eps_c, nG=NG):
    layers, eps_sub = L.stack_cylinder_metasurface(lam, eps_c, eps_c, PERIOD, DIAM, H_CYL, H_FILM, SPACER)
    s = L.solve(lam, PERIOD / 1000.0, layers, eps_sub=eps_sub, nG=nG)
    return s.R, s.T, s.A


def run_conv():
    au = AU["lerer"]
    lines = ["сходимость по числу гармоник nG, формула MG, таблица Лерера; A при 450/550/650/750 нм"]
    for nG in (61, 121, 221, 361):
        t0 = time.time()
        row = []
        for lam in (450, 550, 650, 750):
            eps_c = L.mg(EPS_H, au(lam), L.C_FILL)
            row.append("%.4f" % solve_grating(lam, eps_c, nG)[2])
        lines.append("nG = %3d (%.0f с): %s" % (nG, time.time() - t0, " ".join(row)))
        print(lines[-1])
    io.open(L.RESULTS / "mixing_rules_on_grating_convergence.txt", "w", encoding="utf-8").write("\n".join(lines) + "\n")


def run_main():
    out = io.open(L.RESULTS / "mixing_rules_on_grating.csv", "w", encoding="utf-8", newline="\n")
    out.write("table,formula,lambda_nm,Re_eps_eff,Im_eps_eff,R,T,A\n")
    for tname, au in AU.items():
        for fname in FORMULAS:
            t0 = time.time()
            As = []
            for lam in LAMS:
                eps_c = L.MIXING[fname](EPS_H, au(lam), L.C_FILL, lam)
                R, T, A = solve_grating(lam, eps_c)
                out.write("%s,%s,%d,%.6f,%.6f,%.6f,%.6f,%.6f\n" % (tname, fname, lam, eps_c.real, eps_c.imag, R, T, A))
                As.append(A)
            out.flush()
            print("%-6s %-10s A(600..750 шаг 50): %s  среднее %.3f  (%.0f с)"
                  % (tname, fname, " ".join("%.3f" % As[LAMS.index(l)] for l in (600, 650, 700, 750)),
                     float(np.mean(As)), time.time() - t0))
    out.close()


def run_extra():
    """Первый порядок и Лоренц-Лоренц (обе таблицы), два сорта частиц и ошибка COMPOSITE_3 (J&C)."""
    out = io.open(L.RESULTS / "mixing_rules_on_grating_extra.csv", "w", encoding="utf-8", newline="\n")
    out.write("table,formula,lambda_nm,Re_eps_eff,Im_eps_eff,R,T,A\n")
    for tname, au in AU.items():
        for fname in EXTRA_FORMULAS:
            for lam in LAMS10:
                eps_c = L.MIXING[fname](EPS_H, au(lam), L.C_FILL, lam)
                R, T, A = solve_grating(lam, eps_c)
                out.write("%s,%s,%d,%.6f,%.6f,%.6f,%.6f,%.6f\n" % (tname, fname, lam, eps_c.real, eps_c.imag, R, T, A))
            out.flush()
            print("%s %s готово" % (tname, fname))
    pops = {
        "Au10": lambda lam: [(L.eps_au_jc(lam), 0.10)],
        "Ag10": lambda lam: [(L.eps_ag_jc(lam), 0.10)],
        "Cu10": lambda lam: [(L.eps_cu_jc(lam), 0.10)],
        "Au5+Ag5": lambda lam: [(L.eps_au_jc(lam), 0.05), (L.eps_ag_jc(lam), 0.05)],
        "Au5+Cu5": lambda lam: [(L.eps_au_jc(lam), 0.05), (L.eps_cu_jc(lam), 0.05)],
    }
    for pname, fn in pops.items():
        for lam in LAMS10:
            eps_c = L.mg_multi(EPS_H, fn(lam))
            R, T, A = solve_grating(lam, eps_c)
            out.write("jc,%s,%d,%.6f,%.6f,%.6f,%.6f,%.6f\n" % (pname, lam, eps_c.real, eps_c.imag, R, T, A))
        out.flush()
        print("jc %s готово" % pname)
    if vc is not None:
        for lam in LAMS10:
            eps_p = L.eps_au_jc(lam).conjugate()   # порт vibr_t работает в конвенции Лерера
            eps_bug = vc.composite3_as_written(EPS_H, EPS_H, eps_p, eps_p, 0.10, 0.0).conjugate()
            R, T, A = solve_grating(lam, eps_bug)
            out.write("jc,COMPOSITE_3_bug_Au10,%d,%.6f,%.6f,%.6f,%.6f,%.6f\n" % (lam, eps_bug.real, eps_bug.imag, R, T, A))
        print("jc COMPOSITE_3 с ошибкой готово")
    out.close()


def summarize():
    import csv
    rows = list(csv.DictReader(io.open(L.RESULTS / "mixing_rules_on_grating.csv", encoding="utf-8")))
    lines = []
    for tname in AU:
        lines.append("таблица золота: %s" % tname)
        A = {f: np.array([float(r["A"]) for r in rows if r["table"] == tname and r["formula"] == f]) for f in FORMULAS}
        lam = np.array(LAMS, dtype=float)
        lines.append("  формула      среднее  A(600)  A(650)  A(700)  A(730)  A(750)  min(400-650)")
        for f in FORMULAS:
            a = A[f]
            g = lambda l: a[LAMS.index(l)]
            lines.append("  %-10s  %.4f  %.3f   %.3f   %.3f   %.3f   %.3f   %.3f"
                         % (f, a.mean(), g(600), g(650), g(700), g(730), g(750), a[lam <= 650].min()))
        stack = np.array([A[f] for f in FORMULAS])
        spread = stack.max(axis=0) - stack.min(axis=0)
        i = int(np.argmax(spread))
        lines.append("  наибольший разброс между формулами: %.3f при %d нм (%s)"
                     % (spread[i], LAMS[i], ", ".join("%s %.3f" % (f, A[f][i]) for f in FORMULAS)))
        first = next((LAMS[j] for j in range(len(LAMS)) if spread[j] > 0.02), None)
        lines.append("  разброс превышает 0,02 начиная с %s нм" % first)
        for f in FORMULAS[1:]:
            d = A[f] - A["MG"]
            j = int(np.argmax(np.abs(d)))
            lines.append("  %-10s - MG: наибольшее отличие %+.3f при %d нм" % (f, d[j], LAMS[j]))
    txt = "\n".join(lines) + "\n"
    print(txt)
    io.open(L.RESULTS / "mixing_rules_on_grating_summary.txt", "w", encoding="utf-8").write(txt)


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    t0 = time.time()
    if which in ("conv", "all"):
        run_conv()
    if which in ("main", "all"):
        run_main()
        summarize()
    if which in ("extra", "all"):
        run_extra()
    if which == "summary":
        summarize()
    print("готово за %.0f с" % (time.time() - t0))
