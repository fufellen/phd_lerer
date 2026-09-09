"""Четыре правила смешивания на столбиковом поглотителе (исходная конструкция Лерера): P(lambda) для
Максвелла-Гарнетта, Бруггемана, MG+MLWA (R = 15 нм) и формулы Лерера - тот же расчёт, что на решётке
цилиндров (mixing_rules_on_grating.py), но на конструкции статьи о столбиках; таблица золота Лерера.

Запуск: python mixing_rules_on_pillars.py -> results/mixing_rules_on_pillars.csv, _summary.txt
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

AU = L.AuModel("lerer")
LAMS = list(range(400, 751, 10))
FORMULAS = ["MG", "Bruggeman", "MLWA15", "Lerer"]
EPS_H = complex(L.N_HOST ** 2, 0.0)


def main():
    t0 = time.time()
    out = io.open(L.RESULTS / "mixing_rules_on_pillars.csv", "w", encoding="utf-8", newline="\n")
    out.write("design,formula,lambda_nm,Re_eps_eff,Im_eps_eff,R,T,A\n")
    res = {}
    for design_name, d in (("original", L.ORIGINAL), ("improved", L.IMPROVED)):
        for f in FORMULAS:
            A = []
            for lam in LAMS:
                mix = lambda eh, ep, frac, _f=f, _lam=float(lam): L.MIXING[_f](eh, ep, frac, _lam)
                layers, eps_sub = L.stack_pillar_absorber(float(lam), AU, mixing=mix, **d)
                s = L.solve(float(lam), d["period_nm"] / 1000.0, layers, eps_sub=eps_sub,
                            nG=221 if d["period_nm"] >= 400 else 121)
                ec = L.MIXING[f](EPS_H, AU(float(lam)), 0.10, float(lam))
                out.write("%s,%s,%d,%.6f,%.6f,%.6f,%.6f,%.6f\n" % (design_name, f, lam, ec.real, ec.imag, s.R, s.T, s.A))
                A.append(s.A)
            out.flush()
            res[(design_name, f)] = np.array(A)
            print("%-9s %-10s A(600/650/700/730/750) = %s  среднее %.4f (%.0f с)"
                  % (design_name, f, " ".join("%.3f" % res[(design_name, f)][LAMS.index(l)] for l in (600, 650, 700, 730, 750)),
                     res[(design_name, f)].mean(), time.time() - t0))
    out.close()
    lines = []
    for design_name in ("original", "improved"):
        stack = np.array([res[(design_name, f)] for f in FORMULAS])
        spread = stack.max(axis=0) - stack.min(axis=0)
        i = int(np.argmax(spread))
        lines.append("%s: наибольший разброс %.3f при %d нм (%s); разброс > 0,02 с %s нм"
                     % (design_name, spread[i], LAMS[i], ", ".join("%s %.3f" % (f, res[(design_name, f)][i]) for f in FORMULAS),
                        next((LAMS[j] for j in range(len(LAMS)) if spread[j] > 0.02), None)))
        for f in FORMULAS:
            a = res[(design_name, f)]
            lines.append("  %-10s среднее %.4f, A(750) %.3f, min %.3f при %d нм" % (f, a.mean(), a[-1], a.min(), LAMS[int(np.argmin(a))]))
    txt = "\n".join(lines) + "\n"
    print(txt)
    io.open(L.RESULTS / "mixing_rules_on_pillars_summary.txt", "w", encoding="utf-8").write(txt)


if __name__ == "__main__":
    main()
