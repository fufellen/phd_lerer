"""Дополнение к pillar_shapes.py: высокие конусы («глаз мотылька»), у которых градиент показателя
растянут на 300 нм, а основание почти смыкается с соседями. Проверяется, спасает ли конус большая
высота, раз при 100...150 нм он проигрывает квадратному столбику.

Запуск: python pillar_shapes_extra.py
Выход: results/pillar_shapes_extra.csv, строки дописываются в results/pillar_shapes_summary.csv
"""
from __future__ import annotations

import io
import sys
import time

import numpy as np

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import lerer_grcwa as L  # noqa: E402
import pillar_shapes as PS  # noqa: E402

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

EXTRA = {
    "cone_h300": "конус на основании столбика высотой 300 нм",
    "motheye_h300": "конус «глаз мотылька»: основание 0,95 периода, высота 300 нм",
    "motheye_h500": "конус «глаз мотылька»: основание 0,95 периода, высота 500 нм",
    "frustum_h300": "усечённый конус (верх - половина основания) высотой 300 нм",
}


def extra_layers(variant: str, design: str):
    d = PS.DESIGNS[design]
    a = d["pillar_nm"] / 1000.0
    P = d["period_nm"] / 1000.0
    dummy = 1.0 + 0j
    if variant == "cone_h300":
        return L.cone_layers(0.3, a, 0.1 * a, dummy, 16, "square", name="cone300")
    if variant == "motheye_h300":
        return L.cone_layers(0.3, 0.95 * P, 0.05 * P, dummy, 16, "square", name="moth300")
    if variant == "motheye_h500":
        return L.cone_layers(0.5, 0.95 * P, 0.05 * P, dummy, 24, "square", name="moth500")
    if variant == "frustum_h300":
        return L.cone_layers(0.3, a, a / 2, dummy, 16, "square", name="frustum300")
    raise ValueError(variant)


def main():
    out = io.open(L.RESULTS / "pillar_shapes_extra.csv", "w", encoding="utf-8", newline="\n")
    out.write("design,variant,lambda_nm,R,T,A_layers,A_total\n")
    summary = []
    for design in ("improved", "original"):
        d = PS.DESIGNS[design]
        for variant in EXTRA:
            t0 = time.time()
            pl = extra_layers(variant, design)
            rows = []
            for lam in PS.LAMS8:
                layers, eps_sub = L.stack_pillar_absorber(float(lam), PS.AU, pillar_layers=pl,
                                                          **{k: v for k, v in d.items() if k not in ("pillar_nm", "h_pil_nm")})
                s = L.solve(float(lam), d["period_nm"] / 1000.0, layers, eps_sub=eps_sub, nG=PS.NG[design])
                rows.append((s.R, s.T, s.A, 1.0 - s.R))
                out.write("%s,%s,%d,%.6f,%.6f,%.6f,%.6f\n" % (design, variant, lam, s.R, s.T, s.A, 1.0 - s.R))
            out.flush()
            A = np.array([r[2] for r in rows]); At = np.array([r[3] for r in rows])
            summary.append("%s,%s,%.4f,%.4f,%.4f,%.4f" % (design, variant, A.mean(), A.min(), At.mean(), At.min()))
            print("%-9s %-14s A8: %s  среднее %.4f min %.3f  (%.0f с)"
                  % (design, variant, " ".join("%.3f" % a for a in A), A.mean(), A.min(), time.time() - t0))
    out.close()
    with io.open(L.RESULTS / "pillar_shapes_summary.csv", "a", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(summary) + "\n")


if __name__ == "__main__":
    main()
