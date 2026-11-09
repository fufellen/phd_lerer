# -*- coding: utf-8 -*-
"""Сводка перебора: лучшая конструкция каждого состава и толщины.

  python ir_summary.py sweep|refine
"""
import csv
import sys

import numpy as np

import ir_common as C

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

name = {"sweep": "ir_sweep.csv", "refine": "ir_refine.csv"}[sys.argv[1]]
rows = list(csv.DictReader(open(C.data_dir() / name, encoding="utf-8")))
print("строк", len(rows))
cols = [c for c in rows[0] if c.startswith("A") and c[1:].isdigit()]
lams = np.array([float(c[1:]) for c in cols])
seen = set()
for r in rows:
    k = (r["combo"], r["h_layer"], r["period"], r["ax"], r["h_pillar"], r["nG"])
    assert k not in seen, k
    seen.add(k)
    a = np.array([float(r[c]) for c in cols])
    assert a.min() > -1e-3 and a.max() < 1.001, (k, a.min(), a.max())
groups = {}
for r in rows:
    groups.setdefault((r["combo"], r["h_layer"], r["nG"]), []).append(r)
print("%-12s %5s %4s | %5s %5s %5s | %7s %6s %11s | ср.750-2000 | мин.A 400-2000 | ср. без столбиков" % ("состав", "h, нм", "nG", "d", "a", "h_p", "ср.A", "ш.0.9", "участок"))
for (combo, h, ng), g in groups.items():
    g.sort(key=lambda r: (float(r["width90"]), float(r["mean_A"])), reverse=True)
    r = g[0]
    a = np.array([float(r[c]) for c in cols])
    mask = lams >= 750
    bare = C.tmm_spectrum(combo, float(h), lams)[:, 2].mean()
    print("%-12s %5s %4s | %5s %5s %5s | %7s %6s %11s | %.3f | %.3f | %.3f" % (
        combo, h, ng, r["period"], r["ax"], r["h_pillar"], r["mean_A"], r["width90"],
        "%s-%s" % (r["lo90"], r["hi90"]) if r["lo90"] else "-", a[mask].mean(), a.min(), bare))
