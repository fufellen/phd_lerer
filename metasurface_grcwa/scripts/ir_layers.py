# -*- coding: utf-8 -*-
"""Две конструкции поглотителя для статьи об ИК-диапазоне: смесь частиц в одном слое и раздельные слои (вариант Лерера 09.11.2026).

Конструкция задаётся словарём: period, ax (сторона столбика), hp (высота), pillar=(материалы, доли, n основы),
layers=[(материалы, доли, n основы, толщина нм), ...] сверху вниз, buffer (нм, n = 1,45), подложка n = 1,45.
Состав композита - формула Клаузиуса-Моссотти (полная Максвелла-Гарнетта) для нескольких сортов.

  python ir_layers.py lerer        воспроизведение расчёта Лерера 09.11.2026 (его библиотека материалов)
  python ir_layers.py sweep [п]    перебор обеих конструкций для пар материалов -> ir_layers_sweep.csv
  python ir_layers.py refine [n]   уточнение лучших -> ir_layers_refine.csv
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
import ir_materials as M  # noqa: E402

L = C.L
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")


def comp(spec, lam):
    names, fracs, nh = spec
    return L.mg_multi(complex(nh * nh, 0.0), [(M.MATERIALS[n](float(lam)), f) for n, f in zip(names, fracs)])


def stack(d, lam):
    layers = []
    if d["hp"] > 0:
        layers.append(L.Pillar(d["hp"] / 1000.0, d["ax"] / 1000.0, d["ax"] / 1000.0, comp(d["pillar"], lam), 1.0 + 0j, "square", 0.0, "pillars"))
    for names, fracs, nh, t in d["layers"]:
        layers.append(L.Uniform(t / 1000.0, comp((names, fracs, nh), lam), "comp"))
    if d.get("buffer", 0) > 0:
        layers.append(L.Uniform(d["buffer"] / 1000.0, C.EPS_SUB, "buffer"))
    return layers


def solve(d, lam, ng=81, nx=120):
    x = float(lam) * 1.00001 if abs(float(lam) - d["period"]) < 1e-6 else float(lam)
    r = L.solve(x, d["period"] / 1000.0, stack(d, x), eps_sub=C.EPS_SUB, nG=ng, nx=nx)
    return r.A, r.R, r.T


def lerer():
    """Его структура: период 400, столбики 250 x 250 x 150 с Cu в основе n = 2,5, слой 2000 Cu, слой 2000 (или 3000) нитрида, подложка; буфера нет."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import ir_plot as P
    lams = np.arange(400.0, 1751.0, 25.0)
    cases = [("HfN", "HfN_LL", 2000.0), ("TaN", "TaN_LL", 2000.0), ("TiN, h = 2000 нм", "TiN_LL", 2000.0), ("ZrN", "ZrN_LL", 2000.0), ("TiN, h = 3000 нм", "TiN_LL", 3000.0)]
    plt.rcParams.update({"font.size": 12})
    fig, ax = plt.subplots(2, 1, figsize=(8.5, 7.5), sharex=True)
    for lab, mat, h in cases:
        d = dict(period=400.0, ax=250.0, hp=150.0, pillar=(("Cu_L",), (0.1,), 2.5), layers=[(("Cu_L",), (0.1,), 1.77, 2000.0), ((mat,), (0.1,), 1.77, h)], buffer=0.0)
        res = np.array([solve(d, x, ng=121) for x in lams])
        A, R = res[:, 0], res[:, 1]
        print("%-18s A(860)=%.3f A(1700)=%.3f R(860)=%.1f дБ R(1700)=%.1f дБ T(1700)=%.3f" % (lab, np.interp(860, lams, A), np.interp(1700, lams, A), 10 * np.log10(np.interp(860, lams, R)), 10 * np.log10(np.interp(1700, lams, R)), np.interp(1700, lams, res[:, 2])))
        ax[0].plot(lams, 10 * np.log10(R), lw=1.6, label=lab)
        ax[1].plot(lams, A, lw=1.8, label=lab)
    ax[0].set_ylabel("R, дБ")
    ax[1].set_ylabel("Поглощение A")
    ax[1].set_xlabel("Длина волны, нм")
    for a in ax:
        a.grid(alpha=0.3)
    ax[1].legend(fontsize=10)
    fig.tight_layout()
    out = P.media() / "Двухслойный вариант Лерера — повторение методом Фурье-мод.png"
    fig.savefig(out, dpi=170)
    print(out)


# пары: (верхний сорт - резонанс в видимой области, нижний - в ИК), набор постоянных
PAIRS = {
    "Cu+TiN": ("Cu_JC", "TiN_R17"), "Cu+ZrN": ("Cu_JC", "ZrN_N11"), "ZrN+AZO": ("ZrN_N11", "AZO_S22"), "TiN+AZO": ("TiN_R17", "AZO_S22"), "Cu+AZO": ("Cu_JC", "AZO_S22"),
    "Cu+TiN_LL": ("Cu_L", "TiN_LL"), "Cu+HfN_LL": ("Cu_L", "HfN_LL"), "ZrN+AZO_LL": ("ZrN_LL", "AZO2_LL"), "TiN+AZO_LL": ("TiN_LL", "AZO2_LL"),
}
H = 5000.0
PERIODS = (400.0, 600.0, 800.0)
RATIOS = (0.5, 0.625, 0.75)
HEIGHTS = (150.0, 300.0, 450.0, 600.0)
LAMS_COARSE = np.arange(400.0, 2001.0, 100.0)
LAMS_FINE = np.arange(400.0, 2001.0, 25.0)


def designs():
    """Конструкция A - смесь (по 0,10 каждого сорта в слое H); B10 - слои по 0,10; B20 - слои по 0,20 (столько же частиц, сколько в A).
    Доля верхнего слоя - 1/4 или 1/2 толщины H; столбики из композита верхнего слоя (в A - из смеси)."""
    out = []
    for pair, (m1, m2) in PAIRS.items():
        out.append(("%s|A" % pair, dict(pillar=((m1, m2), (0.1, 0.1), 1.77), layers=[((m1, m2), (0.1, 0.1), 1.77, H)])))
        for f in (0.1, 0.2):
            for q in (0.25, 0.5):
                h1 = H * q
                out.append(("%s|B%02d_%g" % (pair, int(f * 100), q), dict(pillar=((m1,), (f,), 1.77), layers=[((m1,), (f,), 1.77, h1), ((m2,), (f,), 1.77, H - h1)])))
    return out


def _job(job):
    name, d, lams, ng = job
    res = [solve(d, x, ng=ng, nx=160 if d["period"] >= 600 else 120) for x in lams]
    return job, res


def _run(path, jobs, lams, procs):
    head = ["design", "period", "ax", "hp", "nG", "mean_A", "width90", "lo90", "hi90", "min_A"] + ["A%d" % x for x in lams] + ["R%d" % x for x in lams]
    new = not os.path.exists(path)
    t0 = time.time()
    print("заданий", len(jobs), flush=True)
    with open(path, "a", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        if new:
            w.writerow(head)
        with mp.Pool(procs) as pool:
            for k, (job, res) in enumerate(pool.imap_unordered(_job, jobs)):
                name, d, _, ng = job
                A = [x[0] for x in res]
                R = [x[1] for x in res]
                wd, lo, hi = C.band(lams, A, 0.9)
                w.writerow([name, "%g" % d["period"], "%g" % d["ax"], "%g" % d["hp"], ng, "%.5f" % np.mean(A), "%.1f" % wd, "" if lo is None else "%.1f" % lo,
                            "" if hi is None else "%.1f" % hi, "%.5f" % min(A)] + ["%.5f" % a for a in A] + ["%.5f" % r for r in R])
                f.flush()
                if k % 100 == 0:
                    print(k, "из", len(jobs), "%.0f с" % (time.time() - t0), flush=True)
    print("готово за %.0f с" % (time.time() - t0), flush=True)


def _done(path):
    if not os.path.exists(path):
        return set()
    return {(r["design"], r["period"], r["ax"], r["hp"], r["nG"]) for r in csv.DictReader(open(path, encoding="utf-8"))}


def sweep(procs=3):
    path = C.data_dir() / "ir_layers_sweep.csv"
    done = _done(path)
    jobs = []
    for (name, base), p, q, hp in itertools.product(designs(), PERIODS, RATIOS, HEIGHTS):
        d = dict(base, period=p, ax=p * q, hp=hp, buffer=325.0)
        if (name, "%g" % p, "%g" % (p * q), "%g" % hp, "81") not in done:
            jobs.append((name, d, LAMS_COARSE, 81))
    _run(path, jobs, LAMS_COARSE, procs)


def refine(procs=3, ng=121):
    src = C.data_dir() / "ir_layers_sweep.csv"
    path = C.data_dir() / "ir_layers_refine.csv"
    done = _done(path)
    rows = list(csv.DictReader(open(src, encoding="utf-8")))
    best = {}
    for r in rows:
        k = r["design"]
        key = (float(r["width90"]), float(r["mean_A"]))
        if k not in best or key > best[k][0]:
            best[k] = (key, r)
    base = dict(designs())
    jobs = []
    for k, (_, r) in best.items():
        d = dict(base[k], period=float(r["period"]), ax=float(r["ax"]), hp=float(r["hp"]), buffer=325.0)
        if (k, r["period"], r["ax"], r["hp"], str(ng)) not in done:
            jobs.append((k, d, LAMS_FINE, ng))
    _run(path, jobs, LAMS_FINE, procs)


if __name__ == "__main__":
    mp.freeze_support()
    cmd = sys.argv[1]
    if cmd == "lerer":
        lerer()
    elif cmd == "sweep":
        sweep(*(int(x) for x in sys.argv[2:]))
    elif cmd == "refine":
        refine(*(int(x) for x in sys.argv[2:]))
