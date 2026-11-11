# -*- coding: utf-8 -*-
"""Рисунки заметки об ИК-диапазоне. Папка вывода - _media заметки в вауте.

  python ir_plot.py materials   постоянные и ослабление за проход
  python ir_plot.py scheme      схема конструкции
  python ir_plot.py spectra     спектры лучших конструкций (нужен ir_refine.csv)
  python ir_plot.py field COMBO h_layer period ax h_pillar lam1 [lam2 ...]   карта |E|^2 в разрезе
"""
import csv
import os
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from grcwa.fft_funs import get_ifft  # noqa: E402

import ir_common as C  # noqa: E402
import ir_materials as M  # noqa: E402
L = C.L

NOTE = "Статья 2 — ИК-диапазон — материалы и подбор конструкций"


# IR_ARTICLE=1: подписи в терминах статьи и главы 4 (метод ОИДУ, метод фурье-мод, МКЭ), файл с суффиксом « — статья»
ARTICLE = os.environ.get("IR_ARTICLE") == "1"
_TERMS = (("программа VIBR_T", "метод ОИДУ"), ("программы VIBR_T", "метода ОИДУ"), ("VIBR_T", "ОИДУ"), ("Фурье-мод", "фурье-мод"),
          ("COMSOL − ", "МКЭ − "), (" (библиотека Лерера)", ", таблицы программы ОИДУ"))


def term(s):
    if ARTICLE:
        for a, b in _TERMS:
            s = s.replace(a, b)
        if s == "COMSOL":
            s = "МКЭ (COMSOL)"
    return s


def out_name(name):
    return name.replace(".png", " — статья.png") if ARTICLE else name


def media():
    root = next(p for p in Path.home().glob("Мой диск*/Obsidian/Vault"))
    d = root / "PhD/Scientific_study/Scientific_study/paper_analysis/from_Lerer/от Лерера/Статья 2 ИК" / (NOTE + "_media")
    d.mkdir(parents=True, exist_ok=True)
    return d


plt.rcParams.update({"font.size": 12})
COL = {"TiN_R17": "#1f77b4", "TiN_S22": "#1f77b4", "TiN_P84": "#1f77b4", "ZrN_N11": "#2ca02c", "AZO_S22": "#d62728",
       "Al_R98": "#7f7f7f", "Cu_JC": "#ff7f0e"}
STY = {"TiN_S22": "--", "TiN_P84": ":"}


def materials():
    lams = np.arange(400.0, 2001.0, 5.0)
    fig, ax = plt.subplots(1, 3, figsize=(10, 4.4) if ARTICLE else (12, 4.6))
    eh = L.N_HOST ** 2
    for n in ("TiN_R17", "TiN_S22", "ZrN_N11", "AZO_S22", "Al_R98", "Cu_JC"):
        e = np.array([M.MATERIALS[n](float(x)) for x in lams])
        ax[0].plot(lams, e.real, color=COL[n], ls=STY.get(n, "-"), lw=1.8, label=M.LABEL[n])
        ax[1].plot(lams, e.imag, color=COL[n], ls=STY.get(n, "-"), lw=1.8)
        att = []
        for x in lams:
            nn = L.n_of(M.composite([n], [0.1], float(x)))
            att.append(4.343 * 4 * np.pi * nn.imag / x * 2000.0)
        ax[2].plot(lams, att, color=COL[n], ls=STY.get(n, "-"), lw=1.8)
    ax[0].axhline(-2 * eh, color="k", lw=1.0, ls="-.")
    ax[0].text(900, -2 * eh + 1.1, "Re ε = −6,27", fontsize=11)
    ax[0].set_yscale("symlog", linthresh=10)
    ax[0].set_ylabel("Re ε частицы")
    ax[1].set_yscale("log")
    ax[1].set_ylabel("Im ε частицы")
    ax[2].set_yscale("log")
    ax[2].axhline(10.0, color="k", lw=1.0, ls="-.")
    ax[2].text(1020, 12, "10 дБ", fontsize=11)
    ax[2].set_ylabel("Ослабление за проход 2 мкм, дБ" if ARTICLE else "Ослабление за проход слоя 2 мкм, дБ")
    for a in ax:
        a.set_xlabel("Длина волны, нм")
        a.set_xlim(400, 2000)
        a.grid(alpha=0.3, which="both")
    handles, labels = ax[0].get_legend_handles_labels()
    fig.legend(handles, labels, fontsize=11, loc="lower center", ncol=3, frameon=False)
    fig.tight_layout(rect=(0, 0.17, 1, 1))
    out = media() / out_name("ИК материалы — проницаемость и ослабление.png")
    fig.savefig(out, dpi=170)
    print(out)


def scheme():
    fig, (a, b) = plt.subplots(1, 2, figsize=(11.5, 4.2), gridspec_kw=dict(width_ratios=[1.75, 1]))
    fig.suptitle("Конструкция: поперечное сечение (слева) и вид сверху (справа)", fontsize=12)
    p, w, hp, hl, hb = 1.0, 0.65, 0.42, 1.25, 0.55
    x0 = 0.0
    a.add_patch(plt.Rectangle((-1.0, -hl - hb - 0.5), 3.0, 0.5, color="#c9d9ea"))
    a.add_patch(plt.Rectangle((-1.0, -hl - hb), 3.0, hb, color="#dbe7f3"))
    a.add_patch(plt.Rectangle((-1.0, -hl), 3.0, hl, color="#c58b6b"))
    for k in range(-1, 2):
        c = x0 + k * p
        a.add_patch(plt.Rectangle((c - w / 2 + 0.5, 0.0), w, hp, color="#c58b6b"))
    a.text(0.5, -hl / 2, "слой композита: матрица n = 1,77 с наночастицами; толщина $h_c$", ha="center", va="center", fontsize=9.5)
    a.text(0.5, -hl - hb / 2, "буфер n = 1,45; 325 нм", ha="center", va="center", fontsize=9.5)
    a.text(0.5, -hl - hb - 0.25, "подложка n = 1,45", ha="center", va="center", fontsize=9.5)
    a.annotate("", xy=(1.0, 0.0), xytext=(1.0, hp), arrowprops=dict(arrowstyle="<->"))
    a.text(1.02, hp / 2, "$h_p$", fontsize=11, va="center")
    a.annotate("", xy=(0.5 - w / 2, hp + 0.12), xytext=(0.5 + w / 2, hp + 0.12), arrowprops=dict(arrowstyle="<->"))
    a.text(0.5, hp + 0.17, "$a$", fontsize=11, ha="center")
    a.annotate("", xy=(0.5, -0.09), xytext=(1.5, -0.09), arrowprops=dict(arrowstyle="<->", color="w"))
    a.text(1.0, -0.24, "$d$", fontsize=11, ha="center", color="w")
    a.annotate("падающая волна", xy=(-0.5, hp + 0.04), xytext=(-0.5, 1.0), ha="center", arrowprops=dict(arrowstyle="-|>"), fontsize=9.5)
    a.set_xlim(-1.0, 2.0)
    a.set_ylim(-hl - hb - 0.5, 1.15)
    a.set_xlabel("x поперёк, схема без соблюдения масштаба")
    a.set_ylabel("z по нормали")
    a.set_xticks([])
    a.set_yticks([])
    for x in range(0, 3):
        for y in range(0, 3):
            b.add_patch(plt.Rectangle((x + 0.175, y + 0.175), 0.65, 0.65, color="#c58b6b"))
    b.set_xlim(0, 3)
    b.set_ylim(0, 3)
    b.set_aspect("equal")
    b.set_xticks([])
    b.set_yticks([])
    b.set_xlabel("x поперёк: период $d$, сторона столбика $a$")
    b.set_ylabel("y поперёк")
    fig.tight_layout()
    out = media() / "ИК конструкция — схема.png"
    fig.savefig(out, dpi=170)
    print(out)


def _rows(name):
    p = C.data_dir() / name
    return list(csv.DictReader(open(p, encoding="utf-8")))


def spectra():
    rows = _rows("ir_refine.csv")
    best = {}
    for r in rows:
        k = r["combo"], r["h_layer"]
        if k not in best or (float(r["width90"]), float(r["mean_A"])) > (float(best[k]["width90"]), float(best[k]["mean_A"])):
            best[k] = r
    cols = [c for c in rows[0] if c.startswith("A") and c[1:].isdigit()]
    lams = np.array([float(c[1:]) for c in cols])
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 5), sharey=True)
    single = ("TiN", "ZrN", "AZO", "Al", "Cu")
    mix = ("Cu+AZO_05", "TiN+AZO_05", "ZrN+AZO_05", "Cu+TiN_05", "Cu+AZO_10", "TiN+AZO_10", "ZrN+AZO_10")
    palette = plt.get_cmap("tab10").colors
    for a, group, title in ((ax[0], single, "Один сорт частиц, f = 0,1"), (ax[1], mix, "Два сорта: число в подписи — доля каждого сорта")):
        for i, c in enumerate(group):
            for h, ls in (("2000", "-"), ("5000", "--")):
                r = best.get((c, h))
                if r is None:
                    continue
                y = [float(r[k]) for k in cols]
                a.plot(lams, y, ls=ls, lw=1.6, color=palette[i],
                       label="%s, %g мкм" % (c.replace("_05", " 0,05").replace("_10", " 0,10"), float(h) / 1000))
        a.axhline(0.9, color="k", lw=0.9, ls="-.")
        a.set_title(title, fontsize=11)
        a.set_xlabel("Длина волны, нм")
        a.set_xlim(400, 2000)
        a.set_ylim(0, 1)
        a.grid(alpha=0.3)
        a.legend(fontsize=9.5, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.15), frameon=False, columnspacing=0.8, handlelength=1.8)
    ax[0].set_ylabel("Поглощение A")
    fig.tight_layout()
    out = media() / "ИК лучшие конструкции — спектры поглощения.png"
    fig.savefig(out, dpi=170, bbox_inches="tight")
    print(out)


def field(combo, h_layer, period, ax, h_p, lams, ng=121, nx=120):
    """Карта |E|^2 в плоскости x-z через центр столбика; нормировка на падающую волну."""
    fig, axes = plt.subplots(1, len(lams), figsize=(4.0 * len(lams), 4.4))
    if len(lams) == 1:
        axes = [axes]
    top_air = 0.25
    for a, lam in zip(axes, lams):
        layers = C.stack(combo, lam, h_layer, ax, h_p, buffer_nm=325.0)
        obj, _ = L.build(lam, period / 1000.0, layers, eps_sub=C.EPS_SUB, nG=ng, nx=nx)
        obj.RT_Solve(normalize=1)
        # |E0|^2 падающей волны: та же решётка без структуры, воздух сверху и снизу
        ref, _ = L.build(lam, period / 1000.0, [L.Uniform(0.1, 1.0 + 0j, "air")], eps_sub=1.0 + 0j, nG=ng, nx=nx)
        ref.RT_Solve(normalize=1)
        fe0, _ = ref.Solve_FieldFourier(0, 0.5)
        e0sq = float(np.mean(sum(np.abs(get_ifft(nx, nx, c, ref.G)) ** 2 for c in fe0)))
        nz_per_um = 60
        zs, xs_img = [], []
        col = []
        # слой 0 - воздух, далее столбики, композит, буфер, подложка
        names = [("воздух", 0, top_air)] + [(l.name, i + 1, l.t_um) for i, l in enumerate(layers)]
        z_plot = 0.0
        for nm, idx, t in names:
            n = max(6, int(t * nz_per_um))
            for z in np.linspace(0.0, t, n, endpoint=False):
                if idx == 0:
                    zrel = 1.0 - top_air + z
                else:
                    zrel = z
                fe, _fh = obj.Solve_FieldFourier(idx, zrel)
                e2 = sum(np.abs(get_ifft(nx, nx, c, obj.G)) ** 2 for c in fe)
                col.append(e2[:, nx // 2])
                zs.append(z_plot + z)
            z_plot += t
        img = np.array(col) / e0sq
        p_um = period / 1000.0
        z_top = top_air
        extent = [-p_um / 2, p_um / 2, zs[-1] - z_top, -z_top]
        im = a.imshow(img, extent=extent, aspect="auto", cmap="inferno", vmin=0, vmax=np.percentile(img, 99.5))
        plt.colorbar(im, ax=a, label="|E|² / |E₀|²")
        hp_um = h_p / 1000.0
        w = ax / 1000.0
        a.plot([-p_um / 2, -w / 2, -w / 2, w / 2, w / 2, p_um / 2], [hp_um, hp_um, 0, 0, hp_um, hp_um], color="c", lw=1.2)
        z_b = hp_um + h_layer / 1000.0
        a.axhline(z_b, color="c", lw=1.2)
        a.set_title("%g нм" % lam, fontsize=11)
        a.set_xlabel("x поперёк, мкм")
        a.set_ylabel("z вглубь от верха столбика, мкм")
    fig.suptitle("Распределение |E|² в поперечном сечении через центр столбика: %s, слой %g мкм, d = %g нм, a = %g нм, $h_p$ = %g нм"
                 % (combo.replace("_", " "), h_layer / 1000, period, ax, h_p), fontsize=11)
    fig.tight_layout()
    out = media() / ("ИК распределение поля — %s.png" % combo.replace("+", "-"))
    fig.savefig(out, dpi=170)
    print(out)


if __name__ == "__main__":
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    cmd = sys.argv[1]
    if cmd == "materials":
        materials()
    elif cmd == "scheme":
        scheme()
    elif cmd == "spectra":
        spectra()
    elif cmd == "field":
        field(sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6]), [float(x) for x in sys.argv[7:]])
