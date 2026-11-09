# -*- coding: utf-8 -*-
"""Поглотители статей 1 и 2 в COMSOL (четверть ячейки, RunAbsorber13) против метода Фурье-мод и программы VIBR_T.

  python ir_comsol_compare.py <файл absorber13_results.csv>
Пишет сводку на экран и рисунок в _media заметки; ответ метода Фурье-мод кэшируется в ir_fourier_cache.csv.
"""
import csv
import os
import sys
from pathlib import Path

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import ir_common as C  # noqa: E402
import ir_plot as P  # noqa: E402

L = C.L
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

D = C.data_dir()
REPO = Path(__file__).resolve().parents[2]
AU = L.AuModel("lerer")
EH = complex(L.N_HOST ** 2, 0.0)
GEOM = {"ZrN+AZO_10_5000": ("ZrN+AZO_10", 800.0, 500.0, 300.0, 5000.0), "TiN+AZO_10_5000": ("TiN+AZO_10", 800.0, 500.0, 300.0, 5000.0)}


def eps_of(label, lam):
    if label.startswith("Au_MG"):
        return L.mg(EH, AU(lam), 0.10)
    return C.eps_combo(GEOM[label][0], lam)


CACHE = D / "ir_fourier_cache.csv"


def fourier(label, lam, ng=221):
    """A методом Фурье-мод с теми же постоянными и геометрией; кэш по (метка, длина волны, nG)."""
    cache = {}
    if CACHE.exists():
        for r in csv.DictReader(open(CACHE, encoding="utf-8")):
            cache[(r["label"], float(r["lambda_nm"]), int(r["nG"]))] = float(r["A"])
    k = (label, float(lam), ng)
    if k in cache:
        return cache[k]
    if label.startswith("Au_MG"):
        p, a, hp, hc = 400.0, 250.0, 100.0, 2000.0
    else:
        _, p, a, hp, hc = GEOM[label]
    e = eps_of(label, lam)
    layers = [L.Pillar(hp / 1000.0, a / 1000.0, a / 1000.0, e, 1.0 + 0j, "square", 0.0, "pillars"), L.Uniform(hc / 1000.0, e, "comp"), L.Uniform(0.325, C.EPS_SUB, "buf")]
    A = L.solve(float(lam), p / 1000.0, layers, eps_sub=C.EPS_SUB, nG=ng, nx=160).A
    new = not CACHE.exists()
    with open(CACHE, "a", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        if new:
            w.writerow(["label", "lambda_nm", "nG", "A"])
        w.writerow([label, "%g" % lam, ng, "%.6f" % A])
    return A


def vibr():
    v = {}
    p = D / "ir_vibr_spectra.csv"
    if p.exists():
        for r in csv.DictReader(open(p, encoding="utf-8")):
            v.setdefault(r["combo"] + "_5000", {})[float(r["lambda_nm"])] = float(r["A"])
    p2 = REPO / "vibr_t" / "results" / "port_c_rect_spectra.csv"
    if p2.exists():
        for r in csv.DictReader(open(p2, encoding="utf-8")):
            if r["series"] == "MG":
                v.setdefault("Au_MG_2000", {})[float(r["lambda_nm"])] = float(r["A"])
    return v


def main(path):
    last = {}
    for r in csv.DictReader(open(path, encoding="utf-8")):
        last[(r["label"], r["case"], r["lambda_nm"])] = r  # при повторном счёте берётся последняя строка
    rows = list(last.values())
    print("строк", len(rows))
    print("--- контроль: матрица переноса")
    for r in rows:
        if r["case"] not in ("slab", "empty"):
            continue
        lam = float(r["lambda_nm"])
        if r["case"] == "slab":
            n = L.n_of(C.eps_combo("ZrN+AZO_10", lam))
            R, T, A = L.tmm_rta(lam, [(n, float(r["h_comp_nm"])), (L.N_SUB, 325.0)], 1.0, L.N_SUB)
        else:
            R, T, A = L.tmm_rta(lam, [(1.0 + 0j, float(r["h_comp_nm"])), (L.N_SUB, 325.0)], 1.0, L.N_SUB)
        print("%-6s %6.0f нм: COMSOL A=%.4f T=%.4f R=%.4f | матрица A=%.4f T=%.4f R=%.4f | элементов %s" % (r["case"], lam, float(r["A"]), float(r["T"]), float(r["R"]), A, T, R, r["elements"]))
    v = vibr()
    by = {}
    fine = {}
    for r in rows:
        if r["case"] == "pillar" and r["label"].endswith("_fine"):
            fine[(r["label"][:-5], float(r["lambda_nm"]))] = r
        elif r["case"] == "pillar" and not r["label"].endswith("_saved"):
            by.setdefault(r["label"], []).append(r)
    if fine:
        print("--- сходимость по сетке: размеры x 0,7")
        base = {(r["label"], float(r["lambda_nm"])): r for r in rows if r["case"] == "pillar"}
        worst = 0.0
        for (lab, x), rf in sorted(fine.items()):
            rb = base.get((lab, x))
            if rb is None:
                continue
            d = float(rf["A"]) - float(rb["A"])
            worst = max(worst, abs(d))
            print("%-18s %6.0f нм: A %.4f (%s эл.) -> %.4f (%s эл.), разность %+.4f" % (lab, x, float(rb["A"]), rb["elements"], float(rf["A"]), rf["elements"], d))
        print("наибольшая разность %.4f" % worst)
    summary = {}
    for label, g in by.items():
        g.sort(key=lambda r: float(r["lambda_nm"]))
        print("---", label)
        lam, ac, af, av = [], [], [], []
        for r in g:
            x = float(r["lambda_nm"])
            lam.append(x)
            ac.append(float(r["A"]))
            af.append(fourier(label, x))
            av.append(v.get(label, {}).get(x, float("nan")))
            print("%6.0f нм: COMSOL %.4f  Фурье-мод %.4f  VIBR_T %s  | C-Ф %+.4f  C-V %s  T %.4f  элементов %s" % (x, ac[-1], af[-1], "%.4f" % av[-1] if av[-1] == av[-1] else "  -   ", ac[-1] - af[-1], "%+.4f" % (ac[-1] - av[-1]) if av[-1] == av[-1] else "  -   ", float(r["T"]), r["elements"]))
        lam, ac, af, av = map(np.array, (lam, ac, af, av))
        d = ac - af
        print("  макс |COMSOL-Фурье| %.4f при %g нм, ср.кв. %.4f; мин. A COMSOL %.4f при %g нм" % (np.abs(d).max(), lam[np.abs(d).argmax()], float(np.sqrt((d ** 2).mean())), ac.min(), lam[ac.argmin()]))
        vok = ~np.isnan(av)
        if vok.any():
            dv = ac[vok] - av[vok]
            print("  макс |COMSOL-VIBR_T| %.4f при %g нм, точек %d" % (np.abs(dv).max(), lam[vok][np.abs(dv).argmax()], int(vok.sum())))
        summary[label] = (lam, ac, af, av)
    plot(summary)


def plot(s):
    plt.rcParams.update({"font.size": 12})
    n = len(s)
    fig, axes = plt.subplots(2, n, figsize=(4.4 * n, 6.4), squeeze=False, sharex="col")
    for j, (label, (lam, ac, af, av)) in enumerate(s.items()):
        title = {"Au_MG": "Au 0,10; слой 2 мкм", "ZrN+AZO_10": "ZrN + AZO по 0,10; 5 мкм", "TiN+AZO_10": "TiN + AZO по 0,10; 5 мкм"}
        key = next(k for k in title if label.startswith(k))
        a, b = axes[0][j], axes[1][j]
        a.plot(lam, af, "-", color="#1f77b4", lw=1.9, label="метод Фурье-мод, 221 гармоника")
        vok = ~np.isnan(av)
        a.plot(lam[vok], av[vok], "o", color="k", ms=4.0, label="программа VIBR_T")
        a.plot(lam, ac, "s", color="#d62728", ms=4.5, mfc="none", label="COMSOL")
        a.set_title(title[key], fontsize=12)
        a.set_ylabel("Поглощение A")
        a.grid(alpha=0.3)
        if j == 0:
            handles, labels = a.get_legend_handles_labels()
        b.axhline(0, color="k", lw=0.8)
        b.plot(lam, ac - af, "s-", color="#d62728", ms=4.0, mfc="none", label="COMSOL − Фурье-мод")
        b.plot(lam[vok], ac[vok] - av[vok], "o-", color="k", ms=3.5, label="COMSOL − VIBR_T")
        b.set_xlabel("Длина волны, нм")
        b.set_ylabel("Разность A")
        b.grid(alpha=0.3)
        b.legend(fontsize=9.5)
    fig.legend(handles, labels, loc="upper center", ncol=3, fontsize=11, frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = P.media() / "COMSOL против метода Фурье-мод и программы VIBR_T.png"
    fig.savefig(out, dpi=170)
    print(out)


if __name__ == "__main__":
    main(sys.argv[1])
