"""Формы столбиков поглотителя: квадрат, круг, эллипс, усечённый конус, конус, две ступени.

Программа Лерера допускает эллиптические и прямоугольные столбики и конусы («Мои предложения»,
пункт 1); в вауте до сих пор считались только квадратные столбики и две ступени. Здесь сечения и
профили столбиков перебираются открытым решателем grcwa для двух конструкций:
- исходной (период 400 нм, столбик 250 x 250 x 100 нм, композит Au 10 % 2000 нм, прослойка 325 нм,
  подложка n = 1,45), сверенной с COMSOL;
- предложенной (период 200 нм, столбик 125 x 125 x 100 нм, композит 10 % 500 нм над 20 % 1500 нм,
  прослойка 50 нм, золотое зеркало), где остаток непоглощённой мощности - отражение от столбиков.
Конус задаётся лестницей из тонких слоёв; сходимость по числу ступенек проверяется отдельно.
Для конструкций с зеркалом приводятся два поглощения: A_layers = 1 - R - T (в конечных слоях, как в
заметке об улучшении) и A_total = 1 - R (вместе с золотом).

Запуск: python pillar_shapes.py [conv|series|fine|all]
Выход: results/pillar_shapes*.csv|txt
"""
from __future__ import annotations

import io
import math
import sys
import time

import numpy as np

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import lerer_grcwa as L  # noqa: E402

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

AU = L.AuModel("lerer")
LAMS8 = [400, 450, 500, 550, 600, 650, 700, 750]
LAMS_FINE = list(range(400, 751, 10))
NG = {"original": 221, "improved": 121}
DESIGNS = {"original": L.ORIGINAL, "improved": L.IMPROVED}
SLICES = 8


def shape_layers(variant: str, a_nm: float, h_nm: float):
    """Список слоёв столбиков (сверху вниз) для варианта формы; композит подставится в stack_pillar_absorber."""
    a, h = a_nm / 1000.0, h_nm / 1000.0
    dummy = 1.0 + 0j
    if variant == "square":
        return [L.Pillar(h, a, a, dummy, shape="square", name="square")]
    if variant == "circle_area":                     # круг той же площади
        d = 2.0 * a / math.sqrt(math.pi)
        return [L.Pillar(h, d, d, dummy, shape="circle", name="circle")]
    if variant == "circle_side":                     # круг диаметром в сторону квадрата
        return [L.Pillar(h, a, a, dummy, shape="circle", name="circle")]
    if variant in ("ellipse_x", "ellipse_y"):        # эллипс той же площади, оси 1,5 : 1
        ay = a * math.sqrt(4.0 / (1.5 * math.pi))
        ax = 1.5 * ay
        return [L.Pillar(h, ax, ay, dummy, shape="ellipse", name="ellipse")]
    if variant == "frustum":                         # усечённый конус: основание a, верх a/2
        return L.cone_layers(h, a, a / 2, dummy, SLICES, "square", name="frustum")
    if variant == "cone":                            # конус: основание a, верх 0,1 a
        return L.cone_layers(h, a, 0.1 * a, dummy, SLICES, "square", name="cone")
    if variant == "cone_h150":                       # конус в полтора раза выше
        return L.cone_layers(1.5 * h, a, 0.1 * a, dummy, int(1.5 * SLICES), "square", name="cone150")
    if variant == "frustum_h150":
        return L.cone_layers(1.5 * h, a, a / 2, dummy, int(1.5 * SLICES), "square", name="frustum150")
    if variant == "two_step":                        # вторая ступень 0,6 a той же высоты сверху (как в COMSOL)
        return [L.Pillar(h, 0.6 * a, 0.6 * a, dummy, shape="square", name="step2"),
                L.Pillar(h, a, a, dummy, shape="square", name="step1")]
    if variant == "cone_circle":                     # круглый конус той же площади основания
        d = 2.0 * a / math.sqrt(math.pi)
        return L.cone_layers(h, d, 0.1 * d, dummy, SLICES, "circle", name="rcone")
    raise ValueError(variant)


VARIANTS = ["square", "circle_area", "circle_side", "ellipse_x", "ellipse_y", "frustum", "cone", "cone_circle",
            "cone_h150", "frustum_h150", "two_step"]
VARIANT_RU = {
    "square": "квадрат (исходная форма)",
    "circle_area": "круг той же площади",
    "circle_side": "круг диаметром в сторону квадрата",
    "ellipse_x": "эллипс 1,5:1 той же площади, E вдоль большой оси",
    "ellipse_y": "эллипс 1,5:1 той же площади, E вдоль малой оси",
    "frustum": "усечённый конус, верх - половина основания",
    "cone": "конус, верх - десятая часть основания",
    "cone_circle": "круглый конус той же площади основания",
    "cone_h150": "конус в полтора раза выше",
    "frustum_h150": "усечённый конус в полтора раза выше",
    "two_step": "две ступени: вторая 0,6 стороны той же высоты",
}


def solve_variant(design: str, variant: str, lam: float, nG=None, slices=None):
    global SLICES
    old = SLICES
    if slices is not None:
        SLICES = slices
    d = DESIGNS[design]
    pl = shape_layers(variant, d["pillar_nm"], d["h_pil_nm"])
    SLICES = old
    layers, eps_sub = L.stack_pillar_absorber(lam, AU, pillar_layers=pl, **{k: v for k, v in d.items()
                                                                            if k not in ("pillar_nm", "h_pil_nm")})
    pol = "y" if variant == "ellipse_y" else "x"
    s = L.solve(lam, d["period_nm"] / 1000.0, layers, eps_sub=eps_sub, nG=nG or NG[design], pol=pol)
    return s.R, s.T, s.A, 1.0 - s.R


def run_conv():
    lines = []
    for nG in (61, 121, 221):
        t0 = time.time()
        row = ["%d: %.4f" % (lam, solve_variant("improved", "square", lam, nG)[2]) for lam in (400, 550, 750)]
        lines.append("предложенная конструкция, квадрат, гармоник %3d (%.0f с): %s" % (nG, time.time() - t0, " ".join(row)))
        print(lines[-1])
    for sl in (4, 8, 16):
        t0 = time.time()
        row = ["%d: %.4f" % (lam, solve_variant("improved", "cone", lam, slices=sl)[2]) for lam in (400, 550, 750)]
        lines.append("предложенная конструкция, конус, ступенек %2d (%.0f с): %s" % (sl, time.time() - t0, " ".join(row)))
        print(lines[-1])
    io.open(L.RESULTS / "pillar_shapes_convergence.txt", "w", encoding="utf-8").write("\n".join(lines) + "\n")


def run_series():
    out = io.open(L.RESULTS / "pillar_shapes_series.csv", "w", encoding="utf-8", newline="\n")
    out.write("design,variant,lambda_nm,R,T,A_layers,A_total\n")
    summary = ["design,variant,mean8_layers,min8_layers,mean8_total,min8_total"]
    for design in DESIGNS:
        for variant in VARIANTS:
            t0 = time.time()
            rows = []
            for lam in LAMS8:
                R, T, A, At = solve_variant(design, variant, lam)
                rows.append((R, T, A, At))
                out.write("%s,%s,%d,%.6f,%.6f,%.6f,%.6f\n" % (design, variant, lam, R, T, A, At))
            out.flush()
            A = np.array([r[2] for r in rows]); At = np.array([r[3] for r in rows])
            summary.append("%s,%s,%.4f,%.4f,%.4f,%.4f" % (design, variant, A.mean(), A.min(), At.mean(), At.min()))
            print("%-9s %-13s A8: %s  среднее %.4f min %.3f  (%.0f с)"
                  % (design, variant, " ".join("%.3f" % a for a in A), A.mean(), A.min(), time.time() - t0))
    out.close()
    io.open(L.RESULTS / "pillar_shapes_summary.csv", "w", encoding="utf-8", newline="\n").write("\n".join(summary) + "\n")


def run_fine():
    out = io.open(L.RESULTS / "pillar_shapes_fine.csv", "w", encoding="utf-8", newline="\n")
    out.write("design,variant,lambda_nm,R,T,A_layers,A_total\n")
    for design, variants in (("improved", ["square", "circle_area", "frustum", "cone", "two_step", "ellipse_x", "ellipse_y"]),
                             ("original", ["square", "cone", "two_step"])):
        for variant in variants:
            t0 = time.time()
            for lam in LAMS_FINE:
                R, T, A, At = solve_variant(design, variant, lam)
                out.write("%s,%s,%d,%.6f,%.6f,%.6f,%.6f\n" % (design, variant, lam, R, T, A, At))
            out.flush()
            print("%s %s: %.0f с" % (design, variant, time.time() - t0))
    out.close()


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    t0 = time.time()
    if which in ("conv", "all"):
        run_conv()
    if which in ("series", "all"):
        run_series()
    if which in ("fine", "all"):
        run_fine()
    print("готово за %.0f с" % (time.time() - t0))
