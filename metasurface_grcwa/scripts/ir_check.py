# -*- coding: utf-8 -*-
"""Сверка лучших ИК-конструкций: сходимость метода фурье-мод по гармоникам и программа Лерера (vibr_rect). Дозапись в ir_check.csv.

  python ir_check.py COMBO h_layer period ax h_pillar lam1 [lam2 ...]
"""
import csv
import os
import sys
import time
from pathlib import Path

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "vibr_t" / "port_c"))
import ir_common as C  # noqa: E402
import vibr_rect as V  # noqa: E402

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

NGS = ((121, 120), (221, 120), (441, 160))
BASES = ((50, 5, 5), (100, 5, 9))
HEAD = "combo h_layer period ax h_pillar lambda_nm method param R T A seconds".split()


def main(combo, h, p, ax, hp, lams):
    path = C.data_dir() / "ir_check.csv"
    new = not path.exists()
    with open(path, "a", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        if new:
            w.writerow(HEAD)
        base = [combo, "%g" % h, "%g" % p, "%g" % ax, "%g" % hp]
        for lam in lams:
            for ng, nx in NGS:
                t0 = time.time()
                r = C.solve_point(combo, lam, p, h, ax, hp, nG=ng, nx=nx)
                s = time.time() - t0
                w.writerow(base + ["%g" % lam, "fourier", ng, "%.6f" % r.R, "%.6f" % r.T, "%.6f" % r.A, "%.1f" % s])
                f.flush()
                print(combo, lam, "фурье nG=%d R=%.5f T=%.5f A=%.5f, %.0f с" % (ng, r.R, r.T, r.A, s), flush=True)
            e = C.eps_combo(combo, lam)
            g = V.Grating(period=p, ax=ax, ay=ax, h_pillar=hp, h_layer=h)
            for b in BASES:
                t0 = time.time()
                r = V.solve(g, lam, e, e, V.Basis(*b), timeout=3000)
                s = time.time() - t0
                w.writerow(base + ["%g" % lam, "vibr", "%d/%d/%d" % b, "%.6f" % r["R"], "%.6f" % r["T"], "%.6f" % r["A"], "%.1f" % s])
                f.flush()
                print(combo, lam, "программа %s R=%.5f T=%.5f A=%.5f, %.0f с" % (b, r["R"], r["T"], r["A"], s), flush=True)


def spectrum(combo, h, p, ax, hp, workers=3):
    """Весь спектр программой Лерера, базис (50, 5, 5), шаг 25 нм -> ir_vibr_spectra.csv."""
    from concurrent.futures import ThreadPoolExecutor

    import numpy as np
    path = C.data_dir() / "ir_vibr_spectra.csv"
    head = "combo h_layer period ax h_pillar lambda_nm R T A seconds".split()
    done = set()
    if path.exists():
        done = {(r["combo"], r["lambda_nm"]) for r in csv.DictReader(open(path, encoding="utf-8"))}
    lams = [x for x in np.arange(400.0, 2001.0, 25.0) if (combo, "%g" % x) not in done]
    g = V.Grating(period=p, ax=ax, ay=ax, h_pillar=hp, h_layer=h)

    def one(lam):
        e = C.eps_combo(combo, lam)
        x = float(lam) * 1.00001 if abs(float(lam) - p) < 1e-6 else float(lam)
        try:
            return lam, V.solve(g, x, e, e, V.Basis(50, 5, 5), timeout=600)
        except RuntimeError:
            # особая точка программы (при d = 800 нм - ровно 600 нм, lambda/d = 3/4): NaN в выводе; отвести на 0,01 %
            print("NaN при %g нм, отвод на 0,01 %%" % x, flush=True)
            return lam, V.solve(g, x * 1.0001, e, e, V.Basis(50, 5, 5), timeout=600)

    new = not path.exists()
    with open(path, "a", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        if new:
            w.writerow(head)
        with ThreadPoolExecutor(workers) as ex:
            for lam, r in ex.map(one, lams):
                w.writerow([combo, "%g" % h, "%g" % p, "%g" % ax, "%g" % hp, "%g" % lam, "%.6f" % r["R"], "%.6f" % r["T"], "%.6f" % r["A"], "%.1f" % r["seconds"]])
                f.flush()
    print("готово", combo, flush=True)


if __name__ == "__main__":
    if sys.argv[1] == "spectrum":
        spectrum(sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6]))
    else:
        main(sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), [float(x) for x in sys.argv[6:]])
