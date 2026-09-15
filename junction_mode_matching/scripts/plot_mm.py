"""Рисунки к сшиванию мод: структура как решена (карта проницаемости на сетке), спектр мод базиса,
сходимость передачи по числу мод, поля на плоскости стыка.

Запуск: python plot_mm.py <tag> [папка вывода]
"""
import csv
import io
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

HERE = os.path.dirname(os.path.abspath(__file__))


RESULTS = os.path.join(HERE, "mm_results.csv")


def load_results(tag):
    rows = [r for r in csv.DictReader(io.open(RESULTS, encoding="utf-8")) if r["tag"] == tag]
    for r in rows:
        # старые таблицы хранили ненормированные потоки и столбцы C1, C2 с потоком падающей моды
        if "C1" in r:
            inc = float(r["C1"]) if r["dir"] == "12" else float(r["biorth1"])  # в старой таблице поток моды 2 попал в этот столбец
            r["balance"] = str(1.0 - (float(r["P_trans"]) + float(r["P_refl"])) / inc)
    return rows


def convergence(tag, dst):
    rows = load_results(tag)
    if not rows:
        print("нет строк для", tag)
        return
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4))
    for direction, color, name in (("12", "#1f77b4", "диэлектрик → плазмон"), ("21", "#d62728", "плазмон → диэлектрик")):
        for variant, ls in ((1, "-"), (2, "--")):
            sel = sorted([(int(r["N"]), float(r["T_guided"]), float(r["balance"]), float(r["R_guided"])) for r in rows
                          if r["dir"] == direction and int(r["variant"]) == variant])
            if not sel:
                continue
            N, T, bal, R = np.array(sel).T
            axes[0].plot(N, T, ls, marker="o", ms=3, color=color, label="%s, проекция %d" % (name, variant))
            axes[1].plot(N, np.abs(bal), ls, marker="o", ms=3, color=color, label="%s, проекция %d" % (name, variant))
    eta = float(rows[0]["eta_conj"])
    axes[0].axhline(eta, color="gray", ls=":", label="первый порядок, перекрытие по потоку %.4f" % eta)
    # сошедшиеся по сетке значения трёхмерной модели COMSOL (reverse_pairs.csv, сетка 0,7) для той же схемы и состояния
    pairs = os.path.join(HERE, "reverse_pairs.csv")
    if os.path.exists(pairs):
        scheme, state = rows[0]["scheme"], rows[0]["state"]
        best = None
        for p in csv.DictReader(io.open(pairs, encoding="utf-8")):
            if p["scheme"] == scheme and p["state"] == state and (best is None or float(p["mesh"]) < float(best["mesh"])):
                best = p
        if best is not None:
            axes[0].axhline(float(best["T_fwd"]), color="#1f77b4", ls="-.", lw=1.2,
                            label="COMSOL, сетка %s: диэлектрик → плазмон %.4f" % (best["mesh"], float(best["T_fwd"])))
            axes[0].axhline(float(best["T_rev"]), color="#d62728", ls="-.", lw=1.2,
                            label="COMSOL, сетка %s: плазмон → диэлектрик %.4f" % (best["mesh"], float(best["T_rev"])))
    axes[0].set_xlabel("число мод базиса на каждой стороне N")
    axes[0].set_ylabel("доля мощности в направляемой моде T")
    axes[0].set_title("сходимость передачи в направляемую моду", fontsize=10)
    axes[0].grid(alpha=0.3)
    axes[0].legend(fontsize=7)
    axes[1].set_xlabel("число мод базиса на каждой стороне N")
    axes[1].set_ylabel("|1 − P_пр − P_отр|, невязка баланса мощности")
    axes[1].set_yscale("log")
    axes[1].set_title("невязка баланса мощности на стыке", fontsize=10)
    axes[1].grid(alpha=0.3)
    axes[1].legend(fontsize=7)
    fig.suptitle("Сшивание мод на стыке %s: сходимость по числу мод" % tag, fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(dst, "mm_convergence_%s.png" % tag), dpi=150)
    print("mm_convergence_%s.png" % tag)


def spectrum(tag, dst):
    rows = list(csv.DictReader(io.open(os.path.join(HERE, "mm_spectrum_%s.csv" % tag), encoding="utf-8")))
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4), sharey=False)
    for ax, side, name in ((axes[0], "1", "сторона 1: подводящий гребень"), (axes[1], "2", "сторона 2: фазовращатель")):
        r = [x for x in rows if x["side"] == side]
        re = np.array([float(x["re_neff"]) for x in r])
        im = np.array([float(x["im_neff"]) for x in r])
        tm = np.array([float(x["tm_fraction"]) for x in r])
        g = np.array([x["guided"] == "1" for x in r])
        sc = ax.scatter(re, im, c=tm, cmap="coolwarm", vmin=0, vmax=1, s=14)
        ax.scatter(re[g], im[g], s=120, facecolors="none", edgecolors="k", label="направляемая мода стыка")
        ax.axvline(1.444, color="gray", ls=":", label="подложка SiO₂, n = 1,444")
        ax.axvline(1.0, color="gray", ls="--", label="воздух, n = 1")
        ax.set_xlabel("Re n_eff")
        ax.set_ylabel("Im n_eff")
        ax.set_title(name, fontsize=10)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7, loc="upper right")
    cb = fig.colorbar(sc, ax=axes, shrink=0.9)
    cb.set_label("доля |E_y|² (квази-TM характер)")
    fig.suptitle("Базис мод сечений на общей сетке (%s): распространяющиеся при Re n > 0, затухающие вдоль оси Im" % tag, fontsize=10)
    fig.savefig(os.path.join(dst, "mm_spectrum_%s.png" % tag), dpi=150, bbox_inches="tight")
    print("mm_spectrum_%s.png" % tag)


def fields(tag, dst):
    d = np.load(os.path.join(HERE, "mm_fields_%s.npz" % tag))
    x, y = d["x"], d["y"]
    tri = mtri.Triangulation(x, y)
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    panels = [
        (axes[0, 0], np.real(d["eps1"]), "сторона 1: Re ε на сетке (структура как решена)", "viridis"),
        (axes[0, 1], np.real(d["eps2"]), "сторона 2: Re ε на сетке (золото обрезано по шкале)", "viridis"),
        (axes[0, 2], np.abs(d["ey_g2"]), "|E_y| плазмонной моды (сторона 2)", "inferno"),
        (axes[1, 0], np.abs(d["ey_inc"]), "|E_y| падающей моды подвода (сторона 1)", "inferno"),
        (axes[1, 1], np.abs(d["ey_tr"]), "|E_y| прошедшего поля при z = 0⁺ (все моды стороны 2)", "inferno"),
        (axes[1, 2], np.abs(d["ey_ref"]), "|E_y| отражённого поля при z = 0⁻ (все моды стороны 1)", "inferno"),
    ]
    for ax, val, title, cmap in panels:
        v = np.asarray(val, dtype=float)
        vmax = np.percentile(v, 99.5) if "ε" not in title else 12.5
        if "ε" in title:
            v = np.clip(v, -1, 12.5)
            vmin = -1
        else:
            vmin = 0
        tp = ax.tripcolor(tri, v, shading="gouraud", cmap=cmap, vmin=vmin, vmax=vmax)
        fig.colorbar(tp, ax=ax, shrink=0.85)
        ax.set_xlim(0, 0.8)
        ax.set_ylim(-0.45, 0.45)
        ax.set_aspect("equal")
        ax.set_xlabel("x поперёк, мкм (x = 0 — плоскость симметрии)")
        ax.set_ylabel("y по высоте, мкм (y = 0 — верх буфера)")
        ax.set_title(title, fontsize=9)
    fig.suptitle("Поперечное сечение (плоскость стыка z = 0), %s" % tag, fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(dst, "mm_fields_%s.png" % tag), dpi=130)
    print("mm_fields_%s.png" % tag)


def main():
    global RESULTS
    tag = sys.argv[1]
    dst = sys.argv[2] if len(sys.argv) > 2 else HERE
    if len(sys.argv) > 3:
        RESULTS = sys.argv[3]
    convergence(tag, dst)
    spectrum(tag, dst)
    fields(tag, dst)


if __name__ == "__main__":
    main()
