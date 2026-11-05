# -*- coding: utf-8 -*-
"""Причина различия варианта с k = 0 и прямой подачи: признак bu_p ставится в Vvod.c по k из входного файла.

Слой задан одним и тем же композитом; меняется только k в его строке: 0 (bu_p = 0) и 1e-9 (bu_p = 1).
Сам k = 1e-9 на проницаемость не влияет, потому что композит её перезаписывает.
  python check_bu_p.py -> results/port_c_bu_p.csv
"""
import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "scripts"))
from vibr_run import Basis, Grating, input_text, run_text  # noqa: E402
import vibr_materials as VM  # noqa: E402

RES = HERE.parent / "results"


def eps_mg(metal, lam, n_m=1.77, c=0.1):
    e = VM.libmat(metal, lam)
    k = (e - n_m ** 2) / (e + 2 * n_m ** 2)
    return n_m ** 2 * (1 + 2 * c * k) / (1 - c * k)


g = Grating(period=500, ax=300, ay=300, h_pillar=100, h_layer=450, h_buffer=325)
RES.mkdir(exist_ok=True)
with open(RES / "port_c_bu_p.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow("lambda_nm R_k0 T_k0 A_k0 R_k1e-9 T_k1e-9 A_k1e-9 R_direct T_direct A_direct".split())
    for lam in (450.0, 550.0, 600.0, 650.0, 700.0, 750.0):
        text = input_text(g, Basis(50, 10, 3, 3), lam, 1.77 ** 2 + 0j, eps_mg(3, lam + 0.0001))
        row = [lam]
        out = []
        for k in ("0", "1e-9", None):
            s = text
            if k is not None:
                s = "\n".join(("/ composite layer / 0 3 0.1 %.6f 1.77 %s" % (g.h_layer, k)) if x.startswith("/ composite layer") else x
                              for x in text.splitlines()) + "\n\n"
            r = next(iter(run_text(s, 300)[-1]["rows"].values()))
            row += list(r)
            out.append("%.4f" % r[2])
        w.writerow(row)
        print("%g нм: A при k=0 %s, при k=1e-9 %s, напрямую %s" % (lam, *out), flush=True)
