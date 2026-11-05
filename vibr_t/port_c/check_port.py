# -*- coding: utf-8 -*-
"""Проверки сборки VIBR_T.

  python check_port.py reference   воспроизведение reference_outputs/'a (входной файл inputs/A.c, вариант 2)
  python check_port.py direct      та же точка с проницаемостью, поданной напрямую через vibr_run
  python check_port.py converge    сходимость по числам базисных функций и гармоник -> results/port_c_convergence.csv
  python check_port.py methods     слой как композит с k = 0, напрямую и метод Фурье-мод -> results/port_c_methods.csv
  python check_port.py inputs      цилиндры статьи 1: способ задания столбика и слоя против Фурье-мод -> results/port_c_inputs.csv
"""
import csv
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "scripts"))
from read_output import read  # noqa: E402
from vibr_run import Basis, Grating, EXE, solve  # noqa: E402
import vibr_composite as VC  # noqa: E402
import vibr_materials as VM  # noqa: E402

REF = HERE.parent / "reference_outputs" / "'a"
RES = HERE.parent / "results"


def eps_mg(metal, lam, n_m=1.77, c=0.1):
    """Формула (1) статьи 1 в конвенции программы, Im eps < 0."""
    e = VM.libmat(metal, lam)
    k = (e - n_m ** 2) / (e + 2 * n_m ** 2)
    return n_m ** 2 * (1 + 2 * c * k) / (1 - c * k)


def reference():
    work = EXE.parent
    (work / "out_ref.txt").unlink(missing_ok=True)
    (work / "A.c").write_bytes((HERE.parent / "inputs" / "A.c").read_bytes())
    text = (work / "A.c").read_text(encoding="cp1251")
    # только второй вариант: он содержит те же точки, что и эталон; первый и третий здесь не нужны
    parts = text.split("/ lambda /")
    (work / "A_ref.txt").write_text("/ lambda /" + parts[2].replace("200 5 751", "200 250 751"), encoding="cp1251")
    subprocess.run([str(EXE)], cwd=work, input="1\nA_ref.txt\nout_ref.txt\nn\n", capture_output=True, text=True,
                   encoding="cp1251", errors="replace", timeout=600)
    ours = read(str(work / "out_ref.txt"))[0]["rows"]
    bad = 0
    for b in read(str(REF))[1:]:
        for k, v in b["rows"].items():
            d = max(abs(x - y) for x, y in zip(v, ours[k]))
            bad += d > 1e-9
            print(k, "эталон", v, "сборка", ours[k], "ОК" if d <= 1e-9 else "РАЗЛИЧИЕ")
    print("различий:", bad)


def direct():
    g = Grating(period=500, ax=300, ay=300, h_pillar=100, h_layer=450, h_buffer=325)
    ref = {k: v for k, v in read(str(REF))[1]["rows"].items() if k[1] == "s"}
    for (lam, _), v in sorted(ref.items()):
        e = eps_mg(3, lam + 0.0001)
        r = solve(g, lam, 1.77 ** 2 + 0j, e, Basis(50, 10, 3, 3))
        print("%g нм: эталон %s, напрямую %s" % (lam, v, r))


def converge():
    """Цилиндры с размерами статьи 1: сходимость по базису и сверка с методом Фурье-мод при двух числах гармоник."""
    sys.path.insert(0, str(HERE.parent.parent / "metasurface_grcwa" / "scripts"))
    import lerer_grcwa as L

    RES.mkdir(exist_ok=True)
    grid_2000 = [(10, 20, 3, 3), (16, 20, 3, 3), (22, 20, 3, 3), (22, 20, 5, 5), (22, 20, 7, 7), (22, 20, 9, 9)]
    grid_800 = [(22, 20, 5, 5), (35, 20, 5, 5), (50, 20, 5, 5)]
    with open(RES / "port_c_convergence.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow("lambda_nm h_layer N_g M_z M_fi M_r R T A A_fourier_221 A_fourier_441 seconds".split())
        for lam in (450.0, 550.0, 650.0, 700.0, 750.0):
            e = eps_mg(3, lam + 0.0001)
            for h, grid in ((2000.0, grid_2000), (800.0, grid_800)):
                g = Grating(h_layer=h)
                layers = [L.Pillar(g.h_pillar / 1000, g.ax / 1000, g.ay / 1000, e.conjugate(), 1.0 + 0j, "circle", 0.0, "cyl"),
                          L.Uniform(h / 1000, e.conjugate(), "film"), L.Uniform(g.h_buffer / 1000, g.n_buffer ** 2 + 0j, "buf")]
                af = [L.solve(lam, g.period / 1000, layers, eps_sub=g.n_sub ** 2 + 0j, nG=n).A for n in (221, 441)]
                print("%g нм, слой %g: Фурье-мод nG=221 %.4f, nG=441 %.4f" % (lam, h, af[0], af[1]), flush=True)
                for n_g, m_z, m_fi, m_r in grid:
                    t0 = time.time()
                    r = solve(g, lam, e, e, Basis(n_g, m_z, m_fi, m_r), timeout=1500)
                    w.writerow([lam, h, n_g, m_z, m_fi, m_r, r[0], r[1], r[2], round(af[0], 5), round(af[1], 5), round(time.time() - t0, 1)])
                    f.flush()
                    print("  ", n_g, m_z, m_fi, m_r, r, "%.0f с" % (time.time() - t0), flush=True)


def inputs():
    """Цилиндры статьи 1, золото, формула (1): как задан столбик и слой во входном файле. Эталон - метод Фурье-мод."""
    sys.path.insert(0, str(HERE.parent.parent / "metasurface_grcwa" / "scripts"))
    import lerer_grcwa as L
    from vibr_run import input_text, run_text

    RES.mkdir(exist_ok=True)
    g = Grating()
    bs = Basis(22, 20, 5, 5)
    with open(RES / "port_c_inputs.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow("lambda_nm A_direct A_layer_k0 A_pillar_comp3 A_both A_fourier_441".split())
        for lam in (450.0, 500.0, 550.0, 600.0, 650.0, 700.0, 750.0):
            e = eps_mg(3, lam + 0.0001)
            text = input_text(g, bs, lam, e, e)
            layer = "/ composite layer / 0 3 0.1 %.6f 1.77 0" % g.h_layer

            def variant(k0, comp):
                out = []
                for s in text.splitlines():
                    if k0 and s.startswith("/ composite layer"):
                        s = layer
                    if comp and s.startswith("/ liner"):
                        s = "/ liner / 0 3 0.1"
                    if comp and s.startswith("/ n_l k_l"):
                        s = "/ n_l k_l / 1.77 0"
                    out.append(s)
                return next(iter(run_text("\n".join(out) + "\n\n", 300)[-1]["rows"].values()))[2]

            layers = [L.Pillar(0.1, 0.25, 0.25, e.conjugate(), 1.0 + 0j, "circle", 0.0, "cyl"),
                      L.Uniform(2.0, e.conjugate(), "film"), L.Uniform(0.325, 1.45 ** 2 + 0j, "buf")]
            row = [variant(False, False), variant(True, False), variant(False, True), variant(True, True),
                   L.solve(lam, 0.4, layers, eps_sub=1.45 ** 2 + 0j, nG=441).A]
            w.writerow([lam] + [round(x, 5) for x in row])
            f.flush()
            print("%g нм: напрямую %.4f, слой k=0 %.4f, столбик через COMPOSITE_3 %.4f, оба %.4f, Фурье-мод %.4f" % (lam, *row), flush=True)


def methods():
    """Поглощение цилиндров эталонной геометрии: слой задан как композит с k = 0 (эталон), напрямую с k > 0, метод Фурье-мод."""
    sys.path.insert(0, str(HERE.parent.parent / "metasurface_grcwa" / "scripts"))
    import lerer_grcwa as L
    from vibr_run import input_text, run_text

    RES.mkdir(exist_ok=True)
    g = Grating(period=500, ax=300, ay=300, h_pillar=100, h_layer=450, h_buffer=325)
    with open(RES / "port_c_methods.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow("lambda_nm N_g M_z M_fi M_r A_composite_k0 A_direct A_fourier seconds".split())
        for lam in (450.0, 550.0, 600.0, 650.0, 700.0, 750.0):
            e = eps_mg(3, lam + 0.0001)
            layers = [L.Pillar(0.1, 0.3, 0.3, 1.77 ** 2 + 0j, 1.0 + 0j, "circle", 0.0, "cyl"),
                      L.Uniform(0.45, e.conjugate(), "film"), L.Uniform(0.325, 1.45 ** 2 + 0j, "buf")]
            af = L.solve(lam, 0.5, layers, eps_sub=1.45 ** 2 + 0j, nG=221).A
            for bs in (Basis(50, 10, 3, 3), Basis(110, 20, 5, 5)):
                t0 = time.time()
                text = input_text(g, bs, lam, 1.77 ** 2 + 0j, e)
                lines = text.splitlines()
                k0 = [("/ composite layer / 0 3 0.1 %.6f 1.77 0" % g.h_layer) if s.startswith("/ composite layer") else s for s in lines]
                try:
                    a0 = next(iter(run_text("\n".join(k0) + "\n\n", 900)[-1]["rows"].values()))[2]
                    a1 = next(iter(run_text(text, 900)[-1]["rows"].values()))[2]
                except Exception as x:
                    print(lam, bs, "сбой:", str(x)[:120], flush=True)
                    continue
                w.writerow([lam, bs.n_g, bs.m_z, bs.m_fi, bs.m_r, a0, a1, round(af, 5), round(time.time() - t0, 1)])
                f.flush()
                print("%g нм %s: с k=0 %.4f, напрямую %.4f, Фурье-мод %.4f, %.0f с" % (lam, (bs.n_g, bs.m_z, bs.m_fi, bs.m_r), a0, a1, af, time.time() - t0), flush=True)


{"reference": reference, "direct": direct, "converge": converge, "methods": methods, "inputs": inputs}[sys.argv[1]]()
