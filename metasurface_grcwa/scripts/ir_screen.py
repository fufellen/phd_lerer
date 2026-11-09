# -*- coding: utf-8 -*-
"""Оценка составов без столбиков (матрица переноса) и время одного расчёта методом фурье-мод.

  python ir_screen.py tmm     таблица для слоёв 2, 3, 5 и 10 мкм
  python ir_screen.py speed   время и сходимость одной точки
"""
import os
import sys
import time

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")

import numpy as np  # noqa: E402

import ir_common as C  # noqa: E402

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")


def tmm():
    print("Без столбиков. A на 400-2000 нм: среднее; ширина полосы A >= 0.8 (начало-конец); A при 600, 1000, 1500, 2000 нм")
    for h in (2000.0, 3000.0, 5000.0, 10000.0):
        print("--- слой %g мкм" % (h / 1000))
        for name in C.COMBOS:
            r = C.tmm_spectrum(name, h)
            A = r[:, 2]
            w = C.band(C.LAMS, A, 0.8)
            pick = [A[list(C.LAMS).index(x)] for x in (600.0, 1000.0, 1500.0, 2000.0)]
            print("%-12s ср.=%.3f  полоса 0.8: %6.0f нм (%s)  %s" % (
                name, A.mean(), w[0], "-" if w[1] is None else "%.0f-%.0f" % (w[1], w[2]), " ".join("%.3f" % a for a in pick)))


def speed():
    for name, lam in (("TiN", 500.0), ("TiN", 1500.0), ("AZO", 1500.0)):
        for ng, nx in ((49, 100), (81, 120), (121, 120), (221, 120), (441, 160)):
            t0 = time.time()
            r = C.solve_point(name, lam, 500.0, 2000.0, 300.0, 200.0, nG=ng, nx=nx)
            print("%-5s %6.0f нм nG=%4d nx=%3d: R=%.5f A=%.5f, %.2f с" % (name, lam, ng, nx, r.R, r.A, time.time() - t0), flush=True)


if __name__ == "__main__":
    {"tmm": tmm, "speed": speed}[sys.argv[1]]()
