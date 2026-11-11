# -*- coding: utf-8 -*-
"""Выжимка результатов статьи 2 (ИК-диапазон): сводка отражения лучших конструкций и рисунки для статьи и главы 4.

  python ir_digest.py table     отражение R лучших конструкций (ir_refine.csv, ir_layers_refine.csv)
  python ir_digest.py scheme    схема двух конструкций: смесь в одном слое и раздельные слои
  python ir_digest.py pairs     A и R двух конструкций для трёх пар материалов
"""
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


def _load(name):
    rows = list(csv.DictReader(open(C.data_dir() / name, encoding="utf-8")))
    cols = [c for c in rows[0] if c[0] == "A" and c[1:].isdigit()]
    lams = np.array([float(c[1:]) for c in cols])
    return rows, cols, lams


def _r_stats(r, cols, lams):
    R = np.array([float(r["R" + c[1:]]) for c in cols])
    i = int(np.argmax(R))
    return R.mean(), 10 * np.log10(R[i]), lams[i], 10 * np.log10(max(R.mean(), 1e-9))


def _best(rows, key):
    best = {}
    for r in rows:
        k = key(r)
        if k is None:
            continue
        v = (float(r["width90"]), float(r["mean_A"]))
        if k not in best or v > best[k][0]:
            best[k] = (v, r)
    return {k: v[1] for k, v in best.items()}


def table():
    rows, cols, lams = _load("ir_refine.csv")
    best = _best(rows, lambda r: (r["combo"], r["h_layer"]) if r["combo"] in C.COMBOS and "_S22" not in r["combo"] and "_P84" not in r["combo"]
                 and "57" not in r["combo"] and not r["combo"].endswith("_L") and r["combo"] != "Cu_R98" else None)
    print("Один слой композита, лучшие по ширине A >= 0,9 (метод фурье-мод, 121 гармоника)")
    print("%-12s %4s | %6s %6s | %7s %10s %7s" % ("состав", "h_c", "ср.A", "ш.0,9", "ср.R", "макс.R дБ", "при нм"))
    for (c, h), r in sorted(best.items(), key=lambda x: (x[0][1], x[0][0])):
        m, mx, at, _ = _r_stats(r, cols, lams)
        print("%-12s %4s | %.3f %6s | %.4f %10.1f %7.0f" % (c, h, float(r["mean_A"]), r["width90"], m, mx, at))
    rows, cols, lams = _load("ir_layers_refine.csv")
    print("\nДве конструкции, слой 5 мкм")
    print("%-12s %-9s | %6s %6s %6s | %7s %10s %7s" % ("пара", "вариант", "ср.A", "ш.0,9", "мин.A", "ср.R", "макс.R дБ", "при нм"))
    best = _best(rows, lambda r: (r["design"].split("|")[0], "A" if r["design"].endswith("|A") else r["design"].split("|")[1][:3]))
    for (p, k), r in sorted(best.items()):
        m, mx, at, _ = _r_stats(r, cols, lams)
        print("%-12s %-9s | %.3f %6s %.3f | %.4f %10.1f %7.0f  (%s)" % (p, k, float(r["mean_A"]), r["width90"], float(r["min_A"]), m, mx, at, r["design"]))


def scheme():
    fig, axs = plt.subplots(1, 2, figsize=(10, 4.3))
    p, w, hp, hl, hb = 1.0, 0.62, 0.40, 1.6, 0.4
    mix, top, bot = "#c58b6b", "#d9a066", "#8c6d9e"
    for a, kind in zip(axs, ("a", "b")):
        a.add_patch(plt.Rectangle((-0.5, -hl - hb - 0.35), 3.0, 0.35, color="#c9d9ea"))
        a.add_patch(plt.Rectangle((-0.5, -hl - hb), 3.0, hb, color="#dbe7f3"))
        if kind == "a":
            a.add_patch(plt.Rectangle((-0.5, -hl), 3.0, hl, color=mix))
            a.text(1.0, -hl / 2, "смесь частиц двух сортов,\nпо $C$ каждого сорта", ha="center", va="center", fontsize=12.5)
            pc = mix
            a.set_title("а) смесь в одном слое", fontsize=14)
        else:
            h1 = hl * 0.35
            a.add_patch(plt.Rectangle((-0.5, -h1), 3.0, h1, color=top))
            a.add_patch(plt.Rectangle((-0.5, -hl), 3.0, hl - h1, color=bot))
            a.text(1.0, -h1 / 2, "сорт 1 (резонанс в видимой области), $h_1$", ha="center", va="center", fontsize=12.5)
            a.text(1.0, -h1 - (hl - h1) / 2, "сорт 2 (резонанс в ИК-области), $h_c - h_1$", ha="center", va="center", fontsize=12.5, color="w")
            pc = top
            a.set_title("б) раздельные слои", fontsize=14)
        for k in range(3):
            c = k * p + 0.5 - w / 2
            if c + w <= 2.5:
                a.add_patch(plt.Rectangle((c, 0.0), w, hp, color=pc))
        a.text(1.0, -hl - hb / 2, "буфер, n = 1,45, 325 нм", ha="center", va="center", fontsize=12)
        a.text(1.0, -hl - hb - 0.175, "подложка, n = 1,45", ha="center", va="center", fontsize=12)
        a.annotate("", xy=(2.62, 0.0), xytext=(2.62, -hl), arrowprops=dict(arrowstyle="<->"))
        a.text(2.66, -hl / 2, "$h_c$", fontsize=13.5, va="center")
        a.annotate("", xy=(2.62, 0.0), xytext=(2.62, hp), arrowprops=dict(arrowstyle="<->"))
        a.text(2.66, hp / 2, "$h_p$", fontsize=13.5, va="center")
        a.annotate("", xy=(0.5 - w / 2, hp + 0.1), xytext=(0.5 + w / 2, hp + 0.1), arrowprops=dict(arrowstyle="<->"))
        a.text(0.5, hp + 0.15, "$a$", fontsize=13.5, ha="center")
        a.annotate("", xy=(0.5, hp + 0.42), xytext=(1.5, hp + 0.42), arrowprops=dict(arrowstyle="<->"))
        a.text(1.0, hp + 0.47, "$d$", fontsize=13.5, ha="center")
        a.set_xlim(-0.5, 2.95)
        a.set_ylim(-hl - hb - 0.35, hp + 0.75)
        a.set_xticks([])
        a.set_yticks([])
        for s in a.spines.values():
            s.set_visible(False)
    fig.tight_layout()
    out = P.media() / "Две конструкции — схема.png"
    fig.savefig(out, dpi=170)
    print(out)


PAIRS3 = (("ZrN+AZO", "ZrN + AZO, опубликованные постоянные"), ("Cu+TiN_LL", "Cu + TiN, таблицы библиотеки [Лерер]"),
          ("TiN+AZO_LL", "TiN + AZO, таблицы библиотеки [Лерер]"))
NAMES = {"A": "смесь в одном слое, по 0,10", "B10": "раздельные слои, по 0,10", "B20": "раздельные слои, по 0,20"}
COLORS = {"A": "#d62728", "B10": "#1f77b4", "B20": "#2ca02c"}
STYLES = {"A": "-", "B10": "--", "B20": "-."}


def pairs():
    rows, cols, lams = _load("ir_layers_refine.csv")
    best = _best(rows, lambda r: (r["design"].split("|")[0], "A" if r["design"].endswith("|A") else r["design"].split("|")[1][:3]))
    plt.rcParams.update({"font.size": 12})
    fig, ax = plt.subplots(2, 3, figsize=(11, 7.6), sharex=True)
    for j, (p, title) in enumerate(PAIRS3):
        for k in ("A", "B10", "B20"):
            r = best[(p, k)]
            A = np.array([float(r[c]) for c in cols])
            R = np.array([float(r["R" + c[1:]]) for c in cols])
            ax[0, j].plot(lams, A, color=COLORS[k], ls=STYLES[k], lw=1.7, label=NAMES[k])
            ax[1, j].plot(lams, 10 * np.log10(np.maximum(R, 1e-5)), color=COLORS[k], ls=STYLES[k], lw=1.5)
        ax[0, j].axhline(0.9, color="k", lw=0.8, ls=":")
        ax[0, j].set_ylim(0.5, 1.0)
        ax[1, j].set_ylim(-40, 0)
        ax[0, j].set_title(title.replace(", таблицы библиотеки [Лерер]", ",\nтаблицы программы ОИДУ").replace(", опубликованные", ",\nопубликованные"), fontsize=12)
        ax[1, j].set_xlabel("Длина волны, нм")
        for a in ax[:, j]:
            a.set_xlim(400, 2000)
            a.grid(alpha=0.3)
    ax[0, 0].set_ylabel("Поглощение A")
    ax[1, 0].set_ylabel("Отражение R, дБ")
    h, l = ax[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=2, frameon=False, fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    out = P.media() / "Две конструкции — три пары, поглощение и отражение.png"
    fig.savefig(out, dpi=170)
    print(out)


def spectra5():
    """Спектры поглощения лучших конструкций со слоем 5 мкм: один сорт частиц и смеси двух сортов."""
    rows, cols, lams = _load("ir_refine.csv")
    best = _best(rows, lambda r: (r["combo"], r["h_layer"]) if r["h_layer"] == "5000" else None)
    plt.rcParams.update({"font.size": 12})
    fig, ax = plt.subplots(1, 2, figsize=(10, 6.0), sharey=True)
    pal = {"TiN": "#1f77b4", "ZrN": "#2ca02c", "AZO": "#d62728", "Al": "#7f7f7f", "Cu": "#ff7f0e"}
    for c in ("TiN", "ZrN", "AZO", "Al", "Cu"):
        r = best[(c, "5000")]
        ax[0].plot(lams, [float(r[k]) for k in cols], color=pal[c], lw=1.8, label="%s 0,10" % c)
    for c, col in (("Cu+AZO", "#ff7f0e"), ("TiN+AZO", "#1f77b4"), ("ZrN+AZO", "#2ca02c")):
        for f, ls in (("10", "-"), ("05", "--")):
            r = best[(c + "_" + f, "5000")]
            ax[1].plot(lams, [float(r[k]) for k in cols], color=col, ls=ls, lw=1.8,
                       label="%s по 0,%s" % (c.replace("+", " + "), f))
    ax[0].set_title("один сорт частиц", fontsize=12)
    ax[1].set_title("два сорта частиц", fontsize=12)
    for a in ax:
        a.axhline(0.9, color="k", lw=0.8, ls=":")
        a.set_xlim(400, 2000)
        a.set_ylim(0, 1)
        a.set_xlabel("Длина волны, нм")
        a.grid(alpha=0.3)
        a.legend(fontsize=11, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2, frameon=False)
    ax[0].set_ylabel("Поглощение A")
    fig.tight_layout()
    out = P.media() / "ИК лучшие конструкции 5 мкм — спектры поглощения.png"
    fig.savefig(out, dpi=170)
    print(out)


{"table": table, "scheme": scheme, "pairs": pairs, "spectra5": spectra5}[sys.argv[1]]()
