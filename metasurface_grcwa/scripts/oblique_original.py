"""Наклонное падение на исходный столбиковый поглотитель Лерера (статья 1 по докладу SFM-2026).

Конструкция ORIGINAL: период 400 нм, столбики 250 x 250 x 100 нм, композит 2000 нм с долей Au 0,10,
прослойка 325 нм, подложка n = 1,45; матрица n = 1,77; правило Максвелла-Гарнетта; таблица золота
Лерера. Плоскость падения xz параллельна стороне столбика (phi = 0); p - поле в плоскости падения,
s - перпендикулярно ей. 221 гармоника, как в расчёте четырёх правил на нормальном падении.

Запуск: python oblique_original.py [--table lerer|jc]
Выход: results/oblique_original[_jc].csv и _summary.txt
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

TABLE = "jc" if "--table" in sys.argv and sys.argv[sys.argv.index("--table") + 1] == "jc" else "lerer"
AU = L.AuModel(TABLE)
SUF = "_jc" if TABLE == "jc" else ""
LAMS = list(range(400, 751, 10))
THETAS = [int(t) for t in sys.argv[sys.argv.index("--thetas") + 1].split(",")] if "--thetas" in sys.argv else [0, 15, 30, 45, 60]
SUF_TH = ("_th" + "_".join(str(t) for t in THETAS)) if "--thetas" in sys.argv else ""
POLS = {"p": "x", "s": "y"}


def main():
    t0 = time.time()
    path = L.RESULTS / f"oblique_original{SUF}{SUF_TH}.csv"
    out = io.open(path, "w", encoding="utf-8", newline="\n")
    out.write("theta_deg,pol,lambda_nm,R,T,A\n")
    res = {}
    for th in THETAS:
        for pn, pol in POLS.items():
            if th == 0 and pn == "s" and (0, "p") in res:
                res[(th, pn)] = res[(0, "p")]
                for lam, (R, T, A) in zip(LAMS, res[(0, "p")]):
                    out.write("%d,%s,%d,%.6f,%.6f,%.6f\n" % (th, pn, lam, R, T, A))
                continue
            rows = []
            for lam in LAMS:
                layers, eps_sub = L.stack_pillar_absorber(float(lam), AU, **L.ORIGINAL)
                s = L.solve(float(lam), 0.4, layers, eps_sub=eps_sub, nG=221, pol=pol, theta_deg=float(th))
                rows.append((s.R, s.T, s.A))
                out.write("%d,%s,%d,%.6f,%.6f,%.6f\n" % (th, pn, lam, s.R, s.T, s.A))
            out.flush()
            res[(th, pn)] = rows
            A = np.array([r[2] for r in rows])
            print("theta %2d %s: A(400..650) min %.3f, A(700) %.3f, A(750) %.3f, среднее %.4f (%.0f с)"
                  % (th, pn, A[:26].min(), A[LAMS.index(700)], A[-1], A.mean(), time.time() - t0), flush=True)
    out.close()
    lines = [f"таблица золота: {TABLE}; MG, C = 0,10; 221 гармоника; плоскость падения вдоль стороны столбика"]
    i650 = LAMS.index(650)
    for th in THETAS:
        Ap = np.array([r[2] for r in res[(th, "p")]])
        As = np.array([r[2] for r in res[(th, "s")]])
        Au = 0.5 * (Ap + As)
        lines.append("theta %2d: среднее 400-750 p %.4f s %.4f неполяр. %.4f; min 400-650 p %.3f s %.3f; "
                     "A(700) p %.3f s %.3f; A(750) p %.3f s %.3f"
                     % (th, Ap.mean(), As.mean(), Au.mean(), Ap[:i650 + 1].min(), As[:i650 + 1].min(),
                        Ap[LAMS.index(700)], As[LAMS.index(700)], Ap[-1], As[-1]))
    txt = "\n".join(lines) + "\n"
    io.open(L.RESULTS / f"oblique_original{SUF}{SUF_TH}_summary.txt", "w", encoding="utf-8").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
