# -*- coding: utf-8 -*-
"""Спектры лучших конструкций: метод фурье-мод (nG = 121) против программы Лерера (50, 5, 5); разброс по наборам постоянных."""
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

D = C.data_dir()
vibr = {}
for r in csv.DictReader(open(D / "ir_vibr_spectra.csv", encoding="utf-8")):
    vibr.setdefault(r["combo"], {})[float(r["lambda_nm"])] = float(r["A"])
four = {}
for name in ("ir_refine.csv", "ir_sens.csv"):
    rows = list(csv.DictReader(open(D / name, encoding="utf-8")))
    cols = [c for c in rows[0] if c[0] == "A" and c[1:].isdigit()]
    for r in rows:
        if r["h_layer"] == "5000" and r["period"] == "800" and r["ax"] == "500" and r["h_pillar"] == "300":
            four[r["combo"]] = np.array([float(r[c]) for c in cols])
lams = np.array([float(c[1:]) for c in cols])

plt.rcParams.update({"font.size": 12})
fig, ax = plt.subplots(1, 3, figsize=(13, 4.6))
for combo, col in (("ZrN+AZO_10", "#d62728"), ("TiN+AZO_10", "#1f77b4")):
    v = np.array([vibr[combo][x] for x in lams])
    f = four[combo]
    d = v - f
    w = C.band(lams, v, 0.9)
    print("%-11s программа: мин A=%.4f при %g нм, ср.=%.4f, ш.0.9=%.0f (%s-%s); Фурье: мин %.4f при %g; |разность| макс=%.4f при %g, ср.кв=%.4f" % (
        combo, v.min(), lams[v.argmin()], v.mean(), w[0], w[1], w[2], f.min(), lams[f.argmin()], np.abs(d).max(), lams[np.abs(d).argmax()], np.sqrt((d ** 2).mean())))
    lab = combo.replace("_10", " по 0,10")
    ax[0].plot(lams, f, color=col, lw=2.0, label=P.term("%s, метод Фурье-мод" % lab))
    ax[0].plot(lams, v, color="k" if combo.startswith("ZrN") else "#555555", lw=0, marker="o", ms=3.2, label=P.term("%s, программа VIBR_T" % lab))
    ax[1].plot(lams, d, color=col, lw=1.8, label=lab)
ax[0].axhline(0.9, color="k", lw=0.9, ls="-.")
ax[0].set_ylim(0.8, 1.0)
ax[0].set_ylabel("Поглощение A")
ax[0].legend(fontsize=8.5, loc="lower right")
ax[1].axhline(0, color="k", lw=0.8)
ax[1].set_ylabel(P.term("A программы VIBR_T − A метода Фурье-мод"))
ax[1].legend(fontsize=9)
for combo, col, lab, ls in (("ZrN+AZO_10", "#d62728", "ZrN (Naik) + AZO 63 нм", "-"), ("ZrN+AZO57_10", "#d62728", "ZrN (Naik) + AZO 57 нм", "--"),
                            ("TiN+AZO_10", "#1f77b4", "TiN (Reddy) + AZO 63 нм", "-"), ("TiN+AZO57_10", "#1f77b4", "TiN (Reddy) + AZO 57 нм", "--"),
                            ("TiN_S22+AZO_10", "#2ca02c", "TiN (Saha) + AZO 63 нм", "-")):
    ax[2].plot(lams, four[combo], color=col, ls=ls, lw=1.7, label=lab)
    w = C.band(lams, four[combo], 0.9)
    print("%-16s мин %.3f при %g, ср.=%.3f, ш.0.9=%.0f (%s-%s)" % (combo, four[combo].min(), lams[four[combo].argmin()], four[combo].mean(), w[0], w[1], w[2]))
ax[2].axhline(0.9, color="k", lw=0.9, ls="-.")
ax[2].set_ylim(0.8, 1.0)
ax[2].set_ylabel("Поглощение A")
ax[2].legend(fontsize=8, loc="lower right")
for a in ax:
    a.set_xlabel("Длина волны, нм")
    a.set_xlim(400, 2000)
    a.grid(alpha=0.3)
fig.tight_layout()
out = P.media() / P.out_name("ИК сверка с программой Лерера и разброс постоянных.png")
fig.savefig(out, dpi=170)
print(out)

print("--- таблицы программы Лерера против публикаций, слой 5 мкм")
rows = list(csv.DictReader(open(D / "ir_sens.csv", encoding="utf-8")))
ref = list(csv.DictReader(open(D / "ir_refine.csv", encoding="utf-8")))
for lab, geo in (("Cu_L", ("500", "312.5", "450")), ("Al_L", ("500", "250", "450"))):
    r = next(x for x in rows if x["combo"] == lab)
    base = lab[:2]
    q = next(x for x in ref if x["combo"] == base and x["h_layer"] == "5000" and (x["period"], x["ax"], x["h_pillar"]) == geo)
    print(lab, "ср.A", r["mean_A"], "ш.0.9", r["width90"], r["lo90"], r["hi90"], "|", base, "ср.A", q["mean_A"], "ш.0.9", q["width90"], q["lo90"], q["hi90"])
