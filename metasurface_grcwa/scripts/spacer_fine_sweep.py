"""Провал поглощения предложенной конструкции около 740 нм, который восьмиточечная серия не видела:
развёртка по толщине прослойки перед золотым зеркалом на сетке 10 нм по длине волны.

Запуск: python spacer_fine_sweep.py
Выход: results/spacer_fine_sweep.csv, results/spacer_fine_sweep_summary.txt
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
SPACERS = [0, 25, 50, 75, 100, 125, 150]
LOWER = [(0.20, 1500.0), (0.20, 2000.0), (0.30, 1500.0)]   # нижний слой композита: доля золота, толщина


def run_case(spacer_nm, lower):
    d = dict(period_nm=200.0, pillar_nm=125.0, h_pil_nm=100.0, layers_comp=((0.10, 500.0), lower), spacer_nm=float(spacer_nm), mirror=True)
    A = []
    for lam in LAMS:
        layers, eps_sub = L.stack_pillar_absorber(float(lam), AU, **d)
        s = L.solve(float(lam), 0.2, layers, eps_sub=eps_sub, nG=121)
        A.append(s.A)
    return np.array(A)


def main():
    t0 = time.time()
    out = io.open(L.RESULTS / "spacer_fine_sweep.csv", "w", encoding="utf-8", newline="\n")
    out.write("lower_fill,lower_t_nm,spacer_nm,lambda_nm,A_layers\n")
    lines = ["нижний слой (доля, нм) | прослойка, нм | среднее по сетке 10 нм | наименьшее и где | среднее по 8 точкам"]
    for lower in LOWER:
        for sp in SPACERS:
            A = run_case(sp, lower)
            for lam, a in zip(LAMS, A):
                out.write("%.2f,%.0f,%d,%d,%.6f\n" % (lower[0], lower[1], sp, lam, a))
            out.flush()
            i8 = [LAMS.index(l) for l in (400, 450, 500, 550, 600, 650, 700, 750)]
            j = int(np.argmin(A))
            lines.append("%.2f / %.0f | %3d | %.4f | %.4f при %d нм | %.4f" % (lower[0], lower[1], sp, A.mean(), A[j], LAMS[j], A[i8].mean()))
            print(lines[-1], "(%.0f с)" % (time.time() - t0))
    out.close()
    io.open(L.RESULTS / "spacer_fine_sweep_summary.txt", "w", encoding="utf-8").write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
