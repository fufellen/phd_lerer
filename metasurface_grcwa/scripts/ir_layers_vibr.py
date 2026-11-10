# -*- coding: utf-8 -*-
"""Сверка двухслойных конструкций с программой VIBR_T (несколько слоёв под столбиками, vibr_rect.solve_layers).

Первая строка - проверка новой записи входного файла: однослойная конструкция ZrN + AZO должна дать прежнее число."""
import os
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "vibr_t" / "port_c"))
import ir_layers as LY  # noqa: E402
import vibr_rect as V  # noqa: E402

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")


def case(d, lam):
    g = V.Grating(period=d["period"], ax=d["ax"], ay=d["ax"], h_pillar=d["hp"], h_buffer=d.get("buffer", 0.0))
    layers = [(t, LY.comp((n, f, nh), lam)) for n, f, nh, t in d["layers"]]
    x = lam * 1.0001 if abs(lam - 600.0) < 1e-6 and d["period"] == 800.0 else lam
    r = V.solve_layers(g, x, LY.comp(d["pillar"], lam), layers, V.Basis(50, 5, 5), timeout=900)
    a, rr, t = LY.solve(d, lam, ng=221, nx=160)
    return lam, r["A"], a, r["R"]


def main():
    base = dict(LY.designs())
    jobs = []
    jobs.append(("регрессия ZrN+AZO смесь", dict(base["ZrN+AZO|A"], period=800.0, ax=500.0, hp=300.0, buffer=325.0), (700.0, 1350.0)))
    lerer = dict(period=400.0, ax=250.0, hp=150.0, pillar=(("Cu_L",), (0.1,), 2.5), layers=[(("Cu_L",), (0.1,), 1.77, 2000.0), (("TiN_LL",), (0.1,), 1.77, 2000.0)], buffer=0.0)
    jobs.append(("Лерер 09.11, TiN 2000", lerer, (860.0, 1300.0, 1700.0)))
    jobs.append(("Cu+TiN_LL слои по 0,10, верх 1,25", dict(base["Cu+TiN_LL|B10_0.25"], period=800.0, ax=500.0, hp=300.0, buffer=325.0), (500.0, 1200.0, 1800.0)))
    jobs.append(("ZrN+AZO слои по 0,20, верх 2,5", dict(base["ZrN+AZO|B20_0.5"], period=800.0, ax=500.0, hp=300.0, buffer=325.0), (500.0, 1200.0, 1800.0)))
    flat = [(name, d, lam) for name, d, lams in jobs for lam in lams]
    with ThreadPoolExecutor(3) as ex:
        for (name, d, lam), res in zip(flat, ex.map(lambda j: case(j[1], j[2]), flat)):
            print("%-36s %6.0f нм: VIBR_T A=%.4f  Фурье-мод (221) A=%.4f  разность %+.4f  R(VIBR_T, порядок 0)=%.4f" % (name, lam, res[1], res[2], res[1] - res[2], res[3]), flush=True)


main()
