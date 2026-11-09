# -*- coding: utf-8 -*-
"""Перебор геометрии столбиков для составов ИК-статьи методом фурье-мод. Дозапись, продолжается той же командой.

  python ir_sweep.py sweep [потоков]             грубый перебор, nG = 81, шаг 100 нм -> ir_sweep.csv в папке данных
  python ir_sweep.py refine [число лучших] [nG]  лучшие конструкции каждого состава, шаг 25 нм -> ir_refine.csv
"""
import csv
import itertools
import os
import sys
import time

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")

import multiprocessing as mp  # noqa: E402

import numpy as np  # noqa: E402

import ir_common as C  # noqa: E402

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

MAIN = ("TiN", "ZrN", "AZO", "Al", "Cu", "Cu+AZO_05", "TiN+AZO_05", "ZrN+AZO_05", "Cu+TiN_05", "Cu+AZO_10", "TiN+AZO_10", "ZrN+AZO_10")
PERIODS = (300.0, 400.0, 500.0, 600.0, 800.0)
RATIOS = (0.5, 0.625, 0.75, 0.875)
HEIGHTS = (100.0, 200.0, 300.0, 450.0, 600.0)
LAYERS = (2000.0, 5000.0)
LAMS_COARSE = np.arange(400.0, 2001.0, 100.0)
LAMS_FINE = np.arange(400.0, 2001.0, 25.0)


def spectrum(job):
    name, h_layer, period, ax, h_p, lams, ng, nx = job
    A, R = [], []
    for lam in lams:
        # период, равный длине волны: порог открытия порядков, отводится как в Vvod.c
        x = float(lam) * 1.00001 if abs(float(lam) - period) < 1e-6 else float(lam)
        r = C.solve_point(name, x, period, h_layer, ax, h_p, nG=ng, nx=nx)
        A.append(r.A)
        R.append(r.R)
    return job, A, R


def _run(path, jobs, lams, procs):
    head = ["combo", "h_layer", "period", "ax", "h_pillar", "nG", "mean_A", "width90", "lo90", "hi90"] + \
           ["A%d" % x for x in lams] + ["R%d" % x for x in lams]
    new = not path.exists()
    t0 = time.time()
    print("заданий", len(jobs), "потоков", procs, flush=True)
    with open(path, "a", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        if new:
            w.writerow(head)
            f.flush()
        with mp.Pool(procs) as pool:
            for k, (job, A, R) in enumerate(pool.imap_unordered(spectrum, jobs)):
                name, h_layer, period, ax, h_p, _, ng, _ = job
                wd, lo, hi = C.band(lams, A, 0.9)
                w.writerow([name, "%g" % h_layer, "%g" % period, "%g" % ax, "%g" % h_p, ng, "%.5f" % float(np.mean(A)), "%.1f" % wd,
                            "" if lo is None else "%.1f" % lo, "" if hi is None else "%.1f" % hi] +
                           ["%.5f" % a for a in A] + ["%.5f" % r for r in R])
                f.flush()
                if k % 50 == 0:
                    print(k, "из", len(jobs), name, "всего %.0f с" % (time.time() - t0), flush=True)
    print("готово за %.0f с" % (time.time() - t0), flush=True)


def _done(path):
    if not path.exists():
        return set()
    with open(path, encoding="utf-8") as f:
        return {(r["combo"], r["h_layer"], r["period"], r["ax"], r["h_pillar"], r["nG"]) for r in csv.DictReader(f)}


def sweep(procs=3):
    path = C.data_dir() / "ir_sweep.csv"
    done = _done(path)
    jobs = []
    for name, h, p, q, hp in itertools.product(MAIN, LAYERS, PERIODS, RATIOS, HEIGHTS):
        ax = p * q
        if (name, "%g" % h, "%g" % p, "%g" % ax, "%g" % hp, "81") not in done:
            jobs.append((name, h, p, ax, hp, LAMS_COARSE, 81, 120))
    _run(path, jobs, LAMS_COARSE, procs)


def top(n, sweep_csv):
    """Лучшие конструкции каждого состава и толщины: по ширине полосы A >= 0.9 и по среднему A."""
    rows = list(csv.DictReader(open(sweep_csv, encoding="utf-8")))
    groups = {}
    for r in rows:
        groups.setdefault((r["combo"], r["h_layer"]), []).append(r)
    out = []
    for key, g in groups.items():
        g.sort(key=lambda r: (float(r["width90"]), float(r["mean_A"])), reverse=True)
        pick = g[:n]
        best_mean = max(g, key=lambda r: float(r["mean_A"]))
        if best_mean not in pick:
            pick.append(best_mean)
        out += pick
    return out


def refine(n=2, ng=121, procs=3):
    src = C.data_dir() / "ir_sweep.csv"
    path = C.data_dir() / "ir_refine.csv"
    done = _done(path)
    jobs = []
    for r in top(n, src):
        key = (r["combo"], r["h_layer"], r["period"], r["ax"], r["h_pillar"], str(ng))
        if key not in done:
            nx = 160 if float(r["period"]) >= 600 else 120
            jobs.append((r["combo"], float(r["h_layer"]), float(r["period"]), float(r["ax"]), float(r["h_pillar"]), LAMS_FINE, ng, nx))
    _run(path, jobs, LAMS_FINE, procs)


def sens(h, period, ax, h_p, names, ng=121, procs=3):
    """Та же конструкция с другими наборами постоянных -> ir_sens.csv."""
    path = C.data_dir() / "ir_sens.csv"
    done = _done(path)
    nx = 160 if period >= 600 else 120
    jobs = [(n, h, period, ax, h_p, LAMS_FINE, ng, nx) for n in names
            if (n, "%g" % h, "%g" % period, "%g" % ax, "%g" % h_p, str(ng)) not in done]
    _run(path, jobs, LAMS_FINE, procs)


if __name__ == "__main__":
    mp.freeze_support()
    cmd = sys.argv[1]
    if cmd == "sweep":
        sweep(*(int(x) for x in sys.argv[2:]))
    elif cmd == "refine":
        refine(*(int(x) for x in sys.argv[2:]))
    elif cmd == "sens":
        sens(float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), sys.argv[6:])
