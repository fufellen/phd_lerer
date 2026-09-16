"""Прямые стыки для контролей каскада: подводящий гребень -> сердцевина без металла -> сердцевина с полной
плёнкой, а также прямой стык гребень -> полная плёнка (соединение встык sync405) и без-металла -> полная
плёнка (металл начинается скачком). Все три сечения решаются в ОДНОМ прогоне и хранятся вместе с их
матрицами Грама: моды одного сечения из разных прогонов различаются фазами и порядком, поэтому смешивать
перекрёстные матрицы этого файла с базисом taper_basis.py нельзя - контроли считаются только на этом наборе.

Запуск: python taper_basis_direct.py [Nmax] [scale]
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from taper_basis import build_mesh, assign, log, K0, LAM, N_SIO2, pick_guided  # noqa: E402
from mode_matching_junction import Side  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    scale = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
    tag = "A_N%d_s%.2f" % (nmax, scale)
    from skfem import Basis, ElementTriP0
    from femwell.maxwell.waveguide import compute_modes
    t0 = time.time()
    msh = os.path.join(HERE, "taper_mesh_direct_%s.msh" % tag)
    mesh, cells = build_mesh(scale, msh)
    basis0 = Basis(mesh, ElementTriP0())
    log("direct interfaces: сетка %d элементов" % mesh.t.shape[1])
    sigma = K0 ** 2 * 1.6 ** 2
    sides = {}
    store = {}
    for w, lab in ((-1.0, "feed"), (0.0, "w000"), (0.450, "w450")):
        eps = assign(basis0, cells, w)
        modes = compute_modes(basis0, eps, wavelength=LAM, num_modes=nmax, order=1, n_guess=1.6)
        s = Side(list(modes), modes[0].basis.dx, sigma)
        sides[lab] = s
        if w < 0:
            cand = [i for i in range(len(s.n)) if s.n[i].real > N_SIO2 and s.tm[i] > 0.5]
            g = max(cand, key=lambda i: s.n[i].real)
        elif w == 0.0:
            g = pick_guided(s, complex(1.967, 0.0))
        else:
            g = pick_guided(s, complex(1.8884, 0.0245))
        store["n_" + lab] = s.n
        store["U_" + lab] = s.U
        store["C_" + lab] = s.C
        store["tm_" + lab] = s.tm
        store["g_" + lab] = g
        log("  %s: %d мод, направляемая n = %.6f%+.5fi (ранг %d), %.0f с" % (lab, nmax, s.n[g].real, s.n[g].imag, g, time.time() - t0))
    for a, b in (("feed", "w000"), ("w000", "w450"), ("feed", "w450")):
        store["X_%s_%s" % (a, b)] = sides[b].cross(sides[a])
        store["Y_%s_%s" % (a, b)] = sides[a].cross(sides[b])
    np.savez_compressed(os.path.join(HERE, "taper_basis_direct_%s.npz" % tag), labels=np.array(["feed", "w000", "w450"]), **store)
    try:
        os.remove(msh)
    except OSError:
        pass
    log("DONE direct interfaces за %.0f с" % (time.time() - t0))


if __name__ == "__main__":
    main()
