# -*- coding: utf-8 -*-
"""Сводка сравнения двух конструкций (смесь в одном слое и раздельные слои) по ir_layers_refine.csv и рисунки."""
import csv
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import ir_common as C  # noqa: E402
import ir_plot as P  # noqa: E402

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

KIND = {"A": "смесь, по 0,10", "B10_0.25": "слои по 0,10, верхний 1,25 мкм", "B10_0.5": "слои по 0,10, верхний 2,5 мкм",
        "B20_0.25": "слои по 0,20, верхний 1,25 мкм", "B20_0.5": "слои по 0,20, верхний 2,5 мкм"}


def load():
    rows = list(csv.DictReader(open(C.data_dir() / "ir_layers_refine.csv", encoding="utf-8")))
    cols = [c for c in rows[0] if c[0] == "A" and c[1:].isdigit()]
    lams = np.array([float(c[1:]) for c in cols])
    out = {}
    for r in rows:
        pair, kind = r["design"].split("|")
        out[(pair, kind)] = dict(r=r, A=np.array([float(r[c]) for c in cols]), R=np.array([float(r["R" + c[1:]]) for c in cols]))
    return lams, out


def main():
    lams, d = load()
    pairs = list(dict.fromkeys(p for p, _ in d))
    print("%-12s %-32s %5s %5s %4s | %6s %6s %13s %6s | %6s" % ("пара", "конструкция", "d", "a", "h_p", "ср.A", "ш.0,9", "участок", "мин.A", "ср.R"))
    best = {}
    for p in pairs:
        for k in KIND:
            v = d.get((p, k))
            if v is None:
                continue
            r = v["r"]
            print("%-12s %-32s %5s %5s %4s | %.3f %6s %13s %.3f | %.3f" % (p, KIND[k], r["period"], r["ax"], r["hp"], float(r["mean_A"]), r["width90"],
                                                                          "%s-%s" % (r["lo90"], r["hi90"]) if r["lo90"] else "-", float(r["min_A"]), v["R"].mean()))
            kind = "A" if k == "A" else k[:3]
            key = (float(r["width90"]), float(r["mean_A"]))
            if (p, kind) not in best or key > best[(p, kind)][0]:
                best[(p, kind)] = (key, k)
        print()
    plt.rcParams.update({"font.size": 11})
    n = len(pairs)
    cols = 3
    rows = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(14, 3.6 * rows), sharex=True, sharey=True)
    colors = {"A": "#d62728", "B10": "#1f77b4", "B20": "#2ca02c"}
    names = {"A": "смесь в одном слое, по 0,10", "B10": "раздельные слои, по 0,10", "B20": "раздельные слои, по 0,20"}
    for i, p in enumerate(pairs):
        a = axes.flat[i]
        for kind in ("A", "B10", "B20"):
            if (p, kind) not in best:
                continue
            k = best[(p, kind)][1]
            a.plot(lams, d[(p, k)]["A"], color=colors[kind], lw=1.7, label=names[kind])
        a.axhline(0.9, color="k", lw=0.8, ls="-.")
        title = P.term(p.replace("_LL", " (библиотека Лерера)"))
        a.set_title(title, fontsize=11)
        a.set_ylim(0.3, 1.0)
        a.grid(alpha=0.3)
    for a in axes.flat[n:]:
        a.axis("off")
    for a in axes[-1]:
        a.set_xlabel("Длина волны, нм")
    for a in axes[:, 0]:
        a.set_ylabel("Поглощение A")
    h, l = axes.flat[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=3, frameon=False, fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = P.media() / P.out_name("Две конструкции — смесь и раздельные слои — поглощение.png")
    fig.savefig(out, dpi=160)
    print(out)

    fig, axes = plt.subplots(rows, cols, figsize=(14, 3.6 * rows), sharex=True, sharey=True)
    for i, p in enumerate(pairs):
        a = axes.flat[i]
        for kind in ("A", "B10", "B20"):
            if (p, kind) not in best:
                continue
            k = best[(p, kind)][1]
            a.plot(lams, 10 * np.log10(np.maximum(d[(p, k)]["R"], 1e-5)), color=colors[kind], lw=1.5, label=names[kind])
        a.set_title(P.term(p.replace("_LL", " (библиотека Лерера)")), fontsize=11)
        a.set_ylim(-40, 0)
        a.grid(alpha=0.3)
    for a in axes.flat[n:]:
        a.axis("off")
    for a in axes[-1]:
        a.set_xlabel("Длина волны, нм")
    for a in axes[:, 0]:
        a.set_ylabel("Отражение R, дБ")
    fig.legend(h, l, loc="upper center", ncol=3, frameon=False, fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = P.media() / P.out_name("Две конструкции — смесь и раздельные слои — отражение.png")
    fig.savefig(out, dpi=160)
    print(out)


main()
