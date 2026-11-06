# -*- coding: utf-8 -*-
"""Столбики статьи 1 по программе для прямоугольных столбиков (vibr_rect) против метода Фурье-мод.

  python check_rect.py converge [потоков]     сходимость по N_g, M_x = M_y, M_z -> results/port_c_rect_convergence.csv
  python check_rect.py fourier [наибольшее nG]   сходимость метода Фурье-мод по числу гармоник -> results/port_c_rect_fourier.csv
  python check_rect.py spectra N_g M_z M_xy [потоков]   спектры семи рядов статьи -> results/port_c_rect_spectra.csv
  python check_rect.py compare <кривые рисунков> <расчёт Фурье-мод>   сводка -> results/port_c_rect_compare.txt

Все файлы дозаписываются: прерванный прогон можно продолжить той же командой.
"""
import csv
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parent / "scripts"), str(HERE.parent.parent / "metasurface_grcwa" / "scripts")]
import lerer_grcwa as L  # noqa: E402
import vibr_materials as VM  # noqa: E402
import vibr_rect as V  # noqa: E402

RES = HERE.parent / "results"
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

AU = L.AuModel("lerer")
EH = complex(L.N_HOST ** 2, 0.0)
METAL = {"Cu": 1, "Ag": 2, "Au": 3}
SERIES = ("MG", "Bruggeman", "MLWA15", "Lerer", "Cu", "Ag", "Au")


def eps_series(s, lam, f=0.1):
    """Проницаемость композита ряда s при Im eps > 0: те же определения, что в pillar_width_check.py."""
    if s in METAL:
        return L.mg(EH, VM.libmat(METAL[s], lam).conjugate(), f)
    return L.MIXING[s](EH, AU(lam), f, lam)


def done_keys(path, names):
    if not path.exists():
        return set()
    with open(path, encoding="utf-8") as f:
        return {tuple(r[n] for n in names) for r in csv.DictReader(f)}


def converge(workers=4):
    RES.mkdir(exist_ok=True)
    path = RES / "port_c_rect_convergence.csv"
    head = "series lambda_nm N_g M_z M_xy R T A M_z_calc MS seconds".split()
    done = done_keys(path, head[:5])
    grid = [(50, 5, m) for m in (3, 5, 7, 9)] + [(n, 5, 5) for n in (20, 30, 100, 150, 200)] + [(50, m, 5) for m in (3, 8, 12)] + \
           [(100, 5, 7), (150, 5, 7), (100, 5, 9)]
    jobs = [("MG", lam, b) for lam in (450.0, 550.0, 650.0, 700.0, 750.0) for b in grid]
    jobs += [("Ag", lam, b) for lam in (450.0, 550.0, 750.0) for b in ((50, 5, 5), (100, 5, 5), (100, 5, 7), (100, 5, 9))]
    jobs = [j for j in jobs if (j[0], "%g" % j[1], str(j[2][0]), str(j[2][1]), str(j[2][2])) not in done]
    jobs.sort(key=lambda j: j[2][2] ** 3 * j[2][1] ** 3 + j[2][0])
    print("заданий", len(jobs), flush=True)
    new = not path.exists()
    with open(path, "a", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        if new:
            w.writerow(head)

        def one(j):
            s, lam, (n, mz, mxy) = j
            e = eps_series(s, lam)
            r = V.solve(V.Grating(), lam, e, e, V.Basis(n, mz, mxy), timeout=3000)
            return j, r

        t0 = time.time()
        with ThreadPoolExecutor(workers) as ex:
            for j, r in ex.map(one, jobs):
                s, lam, (n, mz, mxy) = j
                w.writerow([s, "%g" % lam, n, mz, mxy, r["R"], r["T"], r["A"], r["M_z_calc"], r["MS"], round(r["seconds"], 1)])
                f.flush()
                print(s, lam, n, mz, mxy, "R=%.5f A=%.4f %.0f с (всего %.0f с)" % (r["R"], r["A"], r["seconds"], time.time() - t0), flush=True)


def fourier(max_ng=1681):
    RES.mkdir(exist_ok=True)
    path = RES / "port_c_rect_fourier.csv"
    head = "series lambda_nm nG R T A seconds".split()
    done = done_keys(path, head[:3])
    new = not path.exists()
    with open(path, "a", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        if new:
            w.writerow(head)
        for s, lams in (("MG", (450.0, 550.0, 650.0, 700.0, 750.0)), ("Ag", (450.0, 480.0, 550.0, 750.0))):
            for lam in lams:
                e = eps_series(s, lam)
                layers = [L.Pillar(0.1, 0.25, 0.25, e, 1.0 + 0j, "square", 0.0, "pillars"), L.Uniform(2.0, e, "comp"),
                          L.Uniform(0.325, L.N_SUB ** 2 + 0j, "spacer")]
                for ng, nx in ((221, 120), (441, 120), (841, 160), (1681, 240)):
                    if (s, "%g" % lam, str(ng)) in done or ng > max_ng:
                        continue
                    t0 = time.time()
                    r = L.solve(lam, 0.4, layers, eps_sub=L.N_SUB ** 2 + 0j, nG=ng, nx=nx)
                    w.writerow([s, "%g" % lam, ng, "%.6f" % r.R, "%.6f" % r.T, "%.6f" % r.A, round(time.time() - t0, 1)])
                    f.flush()
                    print(s, lam, ng, "R=%.5f A=%.5f %.0f с" % (r.R, r.A, time.time() - t0), flush=True)


def spectra(n_g, m_z, m_xy, workers=4):
    RES.mkdir(exist_ok=True)
    path = RES / "port_c_rect_spectra.csv"
    head = "series lambda_nm N_g M_z M_xy R T A seconds".split()
    done = done_keys(path, head[:5])
    lams = [float(x) for x in range(400, 751, 5)]
    # Au и MG - одна и та же модель, но таблицы золота задаются разными функциями, поэтому считаются оба ряда
    jobs = [(s, lam) for lam in lams for s in SERIES if (s, "%g" % lam, str(n_g), str(m_z), str(m_xy)) not in done]
    print("заданий", len(jobs), flush=True)
    new = not path.exists()
    with open(path, "a", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        if new:
            w.writerow(head)

        def one(j):
            e = eps_series(*j)
            # 400 нм = период: порог открытия порядков (+-1, 0) в воздухе. Сам Vvod.c уводит с него множителем 1.00001
            lam = j[1] * 1.00001 if j[1] == V.Grating().period else j[1]
            return j, V.solve(V.Grating(), lam, e, e, V.Basis(n_g, m_z, m_xy), timeout=3000)

        t0 = time.time()
        with ThreadPoolExecutor(workers) as ex:
            for k, (j, r) in enumerate(ex.map(one, jobs)):
                w.writerow([j[0], "%g" % j[1], n_g, m_z, m_xy, r["R"], r["T"], r["A"], round(r["seconds"], 1)])
                f.flush()
                if k % 20 == 0:
                    print(k, j, "A=%.4f, всего %.0f с" % (r["A"], time.time() - t0), flush=True)
    print("готово", flush=True)


def compare(curves_csv, sweep_csv):
    """Спектры программы против кривых рисунков статьи и метода Фурье-мод при 221 гармонике (столбик 250 нм) -> results/port_c_rect_compare.txt."""
    import math

    import numpy as np

    prog = {}
    for r in csv.DictReader(open(RES / "port_c_rect_spectra.csv", encoding="utf-8")):
        prog.setdefault(r["series"], {})[float(r["lambda_nm"])] = (float(r["R"]), float(r["T"]), float(r["A"]))
    four = {}
    for r in csv.DictReader(open(sweep_csv, encoding="utf-8")):
        if float(r["pillar_nm"]) == 250.0:
            four.setdefault(r["series"], {})[float(r["lambda_nm"])] = (float(r["R"]), float(r["T"]), float(r["A"]))
    cur = {}
    for r in csv.DictReader(open(curves_csv, encoding="utf-8")):
        cur.setdefault((r["fig"], r["series"]), []).append((float(r["lambda_nm"]), float(r["value"])))

    def on_curve(fig, s, lams):
        pts = sorted(cur[(fig, s)])
        x, y = np.array([p[0] for p in pts]), np.array([p[1] for p in pts])
        keep = [lam for lam in lams if x[0] <= lam <= x[-1]]
        return keep, np.interp(keep, x, y)

    def stats(d):
        d = np.abs(np.asarray(d))
        return math.sqrt(float((d ** 2).mean())), float(d.max())

    out = ["Программа для прямоугольных столбиков, базис N_g = 50, M_z = 5, M_x = M_y = 5, столбик 250 нм. СКО / наибольшее отличие (при какой длине волны).", ""]
    groups = (("четыре модели, золото", ("MG", "Bruggeman", "MLWA15", "Lerer"), "fig3", "fig4"), ("три материала, формула (1)", ("Cu", "Ag", "Au"), "fig5", "fig6"))
    for title, names, fa, fr in groups:
        out.append("--- " + title)
        worst = [0.0, 0.0, 0.0, 0.0]
        for s in names:
            lams = sorted(prog[s])
            a = np.array([prog[s][x][2] for x in lams])
            f = np.array([four[s][x][2] for x in lams])
            d = a - f
            i = int(np.argmax(np.abs(d)))
            rms, mx = stats(d)
            line = "%-10s A: программа - Фурье-мод 221: %.4f / %.4f (%g нм)" % (s, rms, mx, lams[i])
            worst[0] = max(worst[0], mx)
            keep, c = on_curve(fa, s, lams)
            d2 = np.array([prog[s][x][2] for x in keep]) - c
            j = int(np.argmax(np.abs(d2)))
            rms2, mx2 = stats(d2)
            line += "; программа - кривая рисунка: %.4f / %.4f (%g нм, точек %d)" % (rms2, mx2, keep[j], len(keep))
            worst[1] = max(worst[1], mx2)
            keep, c = on_curve(fr, s, lams)
            rdb = np.array([10 * math.log10(prog[s][x][0]) for x in keep])
            d3 = rdb - c
            rms3, mx3 = stats(d3)
            line += "; R, дБ, программа - кривая: %.2f / %.2f" % (rms3, mx3)
            worst[2] = max(worst[2], mx3)
            keep, c = on_curve(fa, s, lams)
            d4 = np.array([four[s][x][2] for x in keep]) - c
            worst[3] = max(worst[3], float(np.abs(d4).max()))
            line += "; Фурье-мод - кривая: %.4f" % float(np.abs(d4).max())
            out.append(line)
        out.append("наибольшее по группе: программа - Фурье-мод %.4f; программа - рисунок %.4f; R программа - рисунок %.2f дБ; Фурье-мод - рисунок %.4f" % tuple(worst))
        out.append("")
    text = "\n".join(out)
    (RES / "port_c_rect_compare.txt").write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "compare":
        compare(sys.argv[2], sys.argv[3])
    if cmd == "converge":
        converge(*(int(x) for x in sys.argv[2:]))
    elif cmd == "fourier":
        fourier(*(int(x) for x in sys.argv[2:]))
    elif cmd == "spectra":
        spectra(*(int(x) for x in sys.argv[2:]))
