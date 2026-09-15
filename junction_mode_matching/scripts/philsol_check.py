"""Сверка показателей опорной пары вторым открытым решателем мод - philsol (векторные конечные разности).

Сечения те же, что в сшивании мод: подводящий гребень Si 405 x 184 нм на буфере Sb2S3 800 x 120 нм и
фазовращатель Sb2S3 | Au 450 x 10 нм | Si 450 x 180 нм, подложка SiO2, воздух сверху, 1,55 мкм.
Сетка равномерная, поэтому шаг по y задан плёнкой (2,5 нм = четыре ячейки на 10 нм), шаг по x крупнее.
Считается полное сечение по ширине, без стенки симметрии: тип стенки в philsol не задаётся явно.

Запуск: python philsol_check.py [dx_nm] [dy_nm]
"""
import io
import math
import sys
import time

import numpy as np

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

LAM = 1.55
K0 = 2 * math.pi / LAM
N_SIO2, N_SI, N_PCM = 1.444, 3.478, 2.712
N_AU = complex(0.6389, 11.1748)
T_PCM, W_PCM, T_RIDGE, D_AU, W_CORE = 0.120, 0.800, 0.180, 0.010, 0.450
W_FEED, T_FEED = 0.405, 0.184


def index_map(side, x, y):
    n = np.ones((len(x), len(y)), dtype=complex)
    X, Y = np.meshgrid(x, y, indexing="ij")
    n[Y < -T_PCM] = N_SIO2
    buf = (Y >= -T_PCM) & (Y < 0.0) & (np.abs(X) < W_PCM / 2)
    n[buf] = complex(N_PCM, 1e-5)
    if side == 1:
        n[(Y >= 0.0) & (Y < T_FEED) & (np.abs(X) < W_FEED / 2)] = N_SI
    else:
        n[(Y >= 0.0) & (Y < D_AU) & (np.abs(X) < W_CORE / 2)] = N_AU
        n[(Y >= D_AU) & (Y < D_AU + T_RIDGE) & (np.abs(X) < W_CORE / 2)] = N_SI
    return n


def main():
    import philsol as ps
    dx = (float(sys.argv[1]) if len(sys.argv) > 1 else 5.0) * 1e-3
    dy = (float(sys.argv[2]) if len(sys.argv) > 2 else 2.5) * 1e-3
    x = np.arange(-0.70, 0.70 + dx / 2, dx)
    y = np.arange(-0.62, 0.50 + dy / 2, dy)
    out = io.open("philsol_results.csv", "a", encoding="utf-8", newline="\n")
    if out.tell() == 0:
        out.write("side,dx_nm,dy_nm,nx,ny,re_neff,im_neff,seconds\n")
    for side, guess in ((1, 1.89), (2, complex(1.888, 0.024))):
        t0 = time.time()
        n = index_map(side, x, y)
        # philsol хранит показатель по трём осям: изотропная среда - три одинаковых слоя
        n3 = np.stack([n, n, n], axis=-1)
        P, mats = ps.core.eigen_build(K0, n3, dx, dy)
        beta, Ex, Ey = ps.solve.solve(P, K0 * guess, neigs=6)
        neff = beta / K0
        # квази-TM мода: наибольшая доля |Ey|^2; среди них ближайшая к затравке
        best = None
        for k in range(len(neff)):
            ex, ey = Ex[k], Ey[k]
            tm = float(np.sum(np.abs(ey) ** 2) / (np.sum(np.abs(ex) ** 2) + np.sum(np.abs(ey) ** 2)))
            nk = complex(neff[k].real, abs(neff[k].imag))
            print("  side %d cand: n = %.6f%+.5fi  TM %.3f" % (side, nk.real, nk.imag, tm), flush=True)
            if tm > 0.5 and (best is None or abs(nk - guess) < abs(best - guess)):
                best = nk
        dt = time.time() - t0
        print("side %d: n_eff = %.6f%+.5fi (dx %.1f nm, dy %.1f nm, %d x %d, %.0f s)"
              % (side, best.real, best.imag, dx * 1e3, dy * 1e3, len(x), len(y), dt), flush=True)
        out.write("%d,%.2f,%.2f,%d,%d,%.7f,%.5e,%.0f\n" % (side, dx * 1e3, dy * 1e3, len(x), len(y), best.real, best.imag, dt))
        out.flush()
    out.close()


if __name__ == "__main__":
    main()
