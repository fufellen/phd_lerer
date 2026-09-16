"""Подбор формы острия: кусочно-линейный монотонный профиль h(u) с тремя внутренними узлами, минимум потерь
перехода в мере трёхмерной модели (T_fwd * exp(alpha_450 L)) по каскаду местных мод (taper_eme.py).

Запуск: python taper_optimize.py <basis.npz> [N] [длины через запятую]
Выход: taper_optimum.csv и строка MT_TABLE для RunTaper.java.
"""
import io
import math
import os
import sys

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from taper_eme import Basis, evaluate  # noqa: E402

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))


def knots_from_params(p):
    """Шесть параметров -> монотонная таблица (u, h) с концами (0,0) и (1,1)."""
    # u-узлы: доли, накапливаемые через softmax-подобное преобразование
    du = np.exp(p[:3])
    du = du / (du.sum() + 1.0)          # три промежутка плюс остаток до 1
    us = np.cumsum(du)
    dh = np.exp(p[3:])
    dh = dh / (dh.sum() + 1.0)
    hs = np.cumsum(dh)
    tab = [(0.0, 0.0)] + list(zip(us, hs)) + [(1.0, 1.0)]
    return tab


def main():
    path = sys.argv[1]
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 100
    lengths = [float(x) for x in sys.argv[3].split(",")] if len(sys.argv) > 3 else [0.5, 1.0]
    basis = Basis(path, N)
    out = io.open(os.path.join(HERE, "taper_optimum.csv"), "w", encoding="utf-8", newline="\n")
    out.write("L_um,T_fwd_eq,T_rev_eq,loss_fwd_dB,table\n")
    for L in lengths:
        best = None
        for seed in range(3):
            rng = np.random.default_rng(seed)
            p0 = rng.normal(0, 0.7, 6)

            def cost(p):
                tab = knots_from_params(p)
                r = evaluate(basis, "table", L, table=tab)
                return -r["T_fwd_eq"]

            res = minimize(cost, p0, method="Nelder-Mead", options=dict(maxiter=250, xatol=2e-3, fatol=1e-6))
            if best is None or res.fun < best.fun:
                best = res
        tab = knots_from_params(best.x)
        r = evaluate(basis, "table", L, table=tab)
        s = ";".join("%.4f:%.4f" % (u, h) for u, h in tab)
        sys.stdout.flush()
        print("L=%.2f: T_fwd_eq=%.4f (%.3f дБ) T_rev_eq=%.4f  профиль %s" % (L, r["T_fwd_eq"], -10 * math.log10(r["T_fwd_eq"]), r["T_rev_eq"], s))
        out.write("%.3f,%.6f,%.6f,%.4f,%s\n" % (L, r["T_fwd_eq"], r["T_rev_eq"], -10 * math.log10(r["T_fwd_eq"]), s))
    out.close()


if __name__ == "__main__":
    main()
