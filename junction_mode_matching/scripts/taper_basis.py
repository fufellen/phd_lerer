"""Базис местных мод перехода по ширине золота и S-матрицы стыков между соседними ширинами.

Сечение перехода: буфер Sb2S3 800 x 120 нм на SiO2, кремниевый гребень 450 x 190 нм на буфере, в его
подошве плёнка золота 10 нм шириной w (0...450 нм) - вне металла кремний лежит прямо на буфере (как и в
двумерной модели RunTaperSections.java и в трёхмерной модели, где металл задан листом под сердцевиной).
Подводящий волновод - гребень 405 x 184 нм на том же буфере без металла (схема sync405).

Все сечения решаются на ОДНОЙ сетке (ячейки по всем границам всех ширин), поэтому поля соседних сечений
живут на одних квадратурных точках и матрицы сшивания считаются без интерполяции. Для каждой ширины
хранятся несопряжённая и сопряжённая матрицы Грама U, C и показатели; для каждой пары соседних ширин -
перекрёстные матрицы X[n,m] = int e_(i+1),n x h_i,m и Y[m,n] = int e_i,m x h_(i+1),n. Этого достаточно,
чтобы каскадом S-матриц (taper_eme.py) посчитать переход ЛЮБОЙ формы острия: монотонный профиль w(z)
проходит все ширины сетки по порядку.

Запуск: python taper_basis.py [Nmax] [scale]
"""
import io
import math
import os
import sys
import time
from collections import OrderedDict

import numpy as np
from shapely.geometry import box

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mode_matching_junction import Side, K0, LAM, N_SIO2, N_SI, EPS_AU, T_PCM, W_PCM  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
N_PCM_A = 2.712
EPS_PCM = complex(N_PCM_A, 1e-5) ** 2
D_AU = 0.010
T_CORE = 0.190          # гребень над буфером там, где металла нет: 190 нм (Au 10 + Si 180)
W_CORE = 0.450
W_FEED, T_FEED = 0.405, 0.184
WIDTHS_NM = [0, 30, 60, 90, 120, 150, 180, 195, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300,
             320, 340, 360, 390, 420, 450]
XMAX, YSUB, YAIR = 1.0, 0.70, 0.70          # расчётная область A трёхмерной модели
LOG = os.path.join(HERE, "taper_basis_status.txt")


def log(msg):
    print(msg, flush=True)
    with io.open(LOG, "a", encoding="utf-8") as f:
        f.write(msg + "\n")


def build_mesh(scale, msh):
    from femwell.mesh import mesh_from_OrderedDict
    from skfem.io.meshio import from_meshio
    xb = sorted(set([0.0] + [w * 1e-3 / 2 for w in WIDTHS_NM if w > 0] + [W_FEED / 2, W_PCM / 2, 0.60]))
    yb = sorted(set([-0.30, -T_PCM, 0.0, D_AU, T_FEED, T_CORE, 0.40]))
    shapes = OrderedDict()
    res = {}
    cells = []
    for i in range(len(xb) - 1):
        for j in range(len(yb) - 1):
            name = "c%d_%d" % (i, j)
            shapes[name] = box(xb[i], yb[j], xb[i + 1], yb[j + 1])
            xc, yc = 0.5 * (xb[i] + xb[i + 1]), 0.5 * (yb[j] + yb[j + 1])
            h = yb[j + 1] - yb[j]
            cells.append((name, xc, yc))
            if h <= 0.012:
                res[name] = {"resolution": 0.003 * scale, "distance": 0.03}
            elif xc < W_PCM / 2 and -T_PCM - 1e-6 < yc < T_CORE + 1e-6:
                res[name] = {"resolution": 0.018 * scale, "distance": 0.10}
            else:
                res[name] = {"resolution": 0.045 * scale, "distance": 0.25}
    shapes["sub"] = box(0.0, -T_PCM - YSUB, XMAX, -0.30)
    shapes["air"] = box(0.0, -T_PCM - YSUB, XMAX, YAIR)
    res["sub"] = {"resolution": 0.10 * scale, "distance": 0.4}
    res["air"] = {"resolution": 0.10 * scale, "distance": 0.4}
    mesh = from_meshio(mesh_from_OrderedDict(shapes, res, default_resolution_max=0.12 * scale, filename=msh))
    return mesh, cells


def eps_section(w_au, xc, yc):
    """w_au < 0: подводящий гребень 405 x 184 без металла; иначе сердцевина 450 x 190 с плёнкой ширины w_au."""
    if yc < -T_PCM:
        return N_SIO2 ** 2
    if yc < 0.0:
        return EPS_PCM if xc < W_PCM / 2 else 1.0
    if w_au < 0:
        return N_SI ** 2 if (xc < W_FEED / 2 and yc < T_FEED) else 1.0
    if xc < W_CORE / 2 and yc < T_CORE:
        if yc < D_AU and xc < w_au / 2:
            return EPS_AU
        return N_SI ** 2
    return 1.0


def assign(basis0, cells, w_au):
    eps = basis0.zeros(dtype=complex)
    eps[basis0.get_dofs(elements="air")] = 1.0
    eps[basis0.get_dofs(elements="sub")] = N_SIO2 ** 2
    for name, xc, yc in cells:
        eps[basis0.get_dofs(elements=name)] = eps_section(w_au, xc, yc)
    return eps


def pick_guided(side, ref):
    """Направляемая квази-TM мода сечения: ближайшая по комплексному n к опорному значению (продолжение)."""
    n = np.array([complex(z.real, abs(z.imag)) for z in side.n])
    return int(np.argmin(np.abs(n - ref)))


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    scale = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
    tag = "A_N%d_s%.2f" % (nmax, scale)
    log("START taper basis %s %s" % (tag, time.strftime("%Y-%m-%d %H:%M:%S")))
    from skfem import Basis, ElementTriP0
    from femwell.maxwell.waveguide import compute_modes
    t0 = time.time()
    msh = os.path.join(HERE, "taper_mesh_%s.msh" % tag)
    mesh, cells = build_mesh(scale, msh)
    basis0 = Basis(mesh, ElementTriP0())
    log("  сетка: %d элементов" % mesh.t.shape[1])
    n_guess = 1.6
    sigma = K0 ** 2 * n_guess ** 2
    # сечения по порядку прохода: подводящий гребень, затем ширины золота от 0 до 450 нм
    sections = [-1.0] + [w * 1e-3 for w in WIDTHS_NM]
    labels = ["feed"] + ["w%03d" % w for w in WIDTHS_NM]
    store = {}
    prev = None
    ref = None
    guided_track = []
    for k, (w, lab) in enumerate(zip(sections, labels)):
        t1 = time.time()
        eps = assign(basis0, cells, w)
        modes = compute_modes(basis0, eps, wavelength=LAM, num_modes=nmax, order=1, n_guess=n_guess)
        side = Side(list(modes), modes[0].basis.dx, sigma)
        store["n_" + lab] = side.n
        store["tm_" + lab] = side.tm
        store["U_" + lab] = side.U
        store["C_" + lab] = side.C
        # направляемая мода: у подводящего гребня - верхняя квази-TM без потерь; далее - продолжение по n
        if w < 0:
            cand = [i for i in range(len(side.n)) if side.n[i].real > N_SIO2 and side.tm[i] > 0.5]
            g = max(cand, key=lambda i: side.n[i].real)
            ref = complex(1.906, 0.0)     # гребень 450 x 190 без металла (femwell, 03.09.2026: 1.906182)
        else:
            g = pick_guided(side, ref)
            ref = complex(side.n[g].real, abs(side.n[g].imag))
        store["g_" + lab] = g
        guided_track.append((lab, side.n[g], side.tm[g]))
        log("  %s: %d мод за %.0f с; направляемая: n = %.6f%+.5fi (TM %.2f, ранг %d), биорт. %.2e, ветвь у %d"
            % (lab, nmax, time.time() - t1, side.n[g].real, side.n[g].imag, side.tm[g], g, side.biorth, int(side.flipped.sum())))
        if prev is not None:
            plab, pside = prev
            store["X_%s_%s" % (plab, lab)] = side.cross(pside)    # int e_(new) x h_(prev)
            store["Y_%s_%s" % (plab, lab)] = pside.cross(side)    # int e_(prev) x h_(new)
        prev = (lab, side)
        if k % 4 == 0:
            np.savez_compressed(os.path.join(HERE, "taper_basis_%s.npz" % tag), labels=np.array(labels[:k + 1]), **store)
    np.savez_compressed(os.path.join(HERE, "taper_basis_%s.npz" % tag), labels=np.array(labels),
                        widths_nm=np.array([-1] + WIDTHS_NM), elements=mesh.t.shape[1], **store)
    with io.open(os.path.join(HERE, "taper_branch_%s.csv" % tag), "w", encoding="utf-8", newline="\n") as f:
        f.write("section,w_au_nm,re_neff,im_neff,alpha_1um,tm_fraction\n")
        for (lab, n, tm), w in zip(guided_track, [-1] + WIDTHS_NM):
            f.write("%s,%d,%.7f,%.5e,%.5f,%.3f\n" % (lab, w, n.real, abs(n.imag), 2 * K0 * abs(n.imag), tm))
    try:
        os.remove(msh)
    except OSError:
        pass
    log("DONE taper basis %s за %.0f с" % (tag, time.time() - t0))


if __name__ == "__main__":
    main()
