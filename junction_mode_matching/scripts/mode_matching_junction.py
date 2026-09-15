"""Сшивание мод на резком стыке кремниевого гребня и плазмонного фазовращателя Sb2S3/Au.

Полуаналитический расчёт стыка: полный набор собственных мод обоих сечений на ОДНОЙ сетке
(femwell, метод конечных элементов, элементы Неделека первого порядка), непрерывность тангенциальных
E и H на плоскости стыка, спроецированная на моды, - двусторонний метод разложения по собственным
модам с одной границей. Ничего не подгоняется: на входе моды, на выходе амплитуды всех мод по обе
стороны стыка, отражение, передача в направляемую моду, доля излучения и невязка баланса мощности.

Считается в обе стороны: диэлектрик -> плазмон (вход в фазовращатель) и плазмон -> диэлектрик
(выход из него). Нормировка несопряжённая (биортогональная), поэтому потери мод учтены честно.

Обозначения. Сторона a - откуда падает мода, сторона b - куда она проходит. Для мод сторон
U_a[m,k] = int (e_am x h_ak).z dS, X[n,m] = int (e_bn x h_am).z dS, Y[m,n] = int (e_am x h_bn).z dS.
Непрерывность E, спроецированная на h_a, и непрерывность H, спроецированная на e_b, дают
    (U_b + X U_a^{-T} X^T) t = 2 X d,   r = U_a^{-T} X^T t - d,
где d - единичный вектор падающей моды, t - амплитуды прошедших мод, r - отражённых.
Второй вариант проекций (E на h_b, H на e_a) даёт
    (U_b^T + Y^T U_a^{-1} Y) t = 2 Y^T d,   r = d - U_a^{-1} Y t;
разница двух вариантов - мера усечения базиса.

Запуск: python mode_matching_junction.py <scheme> <state> <box> <Nmax>
    scheme: abrupt | shared450 | sync405;  state: a | c;  box: A | B;  Nmax: число мод на сторону.
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

HERE = os.path.dirname(os.path.abspath(__file__))
LAM = 1.55
K0 = 2 * math.pi / LAM
N_SIO2, N_SI = 1.444, 3.478
EPS_AU = complex(0.6389, 11.1748) ** 2
N_PCM = {"a": 2.712, "c": 3.308}
T_PCM, W_PCM, T_RIDGE, D_AU, W_CORE = 0.120, 0.800, 0.180, 0.010, 0.450
SCHEMES = {
    "abrupt": dict(shared=False, w_feed=0.736, t_feed=0.240),
    "shared450": dict(shared=True, w_feed=0.450, t_feed=0.184),
    "sync405": dict(shared=True, w_feed=0.405, t_feed=0.184),
}
# окно A - то же, что у трёхмерной модели COMSOL (половина: x в [0, 1], y в [-0.82, 0.70]);
# окно B шире, чтобы видеть, зависит ли ответ от дискретизации континуума
BOXES = {"A": dict(xmax=1.0, ysub=0.70, yair=0.70), "B": dict(xmax=1.5, ysub=1.10, yair=1.00)}
REF_SHIFTER = {"a": complex(1.8874, 0.0241), "c": complex(2.2772, 0.0150)}


def log(msg, fname):
    print(msg, flush=True)
    with io.open(fname, "a", encoding="utf-8") as f:
        f.write(msg + "\n")


def build_mesh(scheme, boxname, scale, msh):
    from femwell.mesh import mesh_from_OrderedDict
    from skfem.io.meshio import from_meshio
    g = SCHEMES[scheme]
    bx = BOXES[boxname]
    y_lo, y_hi = -T_PCM - bx["ysub"], bx["yair"]
    # разбиения: все границы обоих сечений
    xb = sorted(set([0.0, g["w_feed"] / 2, W_CORE / 2, W_PCM / 2, 0.60]))
    y_feed_top = (-T_PCM + g["t_feed"]) if not g["shared"] else g["t_feed"]
    yb = sorted(set([-0.30, -T_PCM, 0.0, D_AU, y_feed_top, D_AU + T_RIDGE, 0.40]))
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
            elif xc < W_PCM / 2 and -T_PCM - 1e-6 < yc < D_AU + T_RIDGE + 1e-6:
                res[name] = {"resolution": 0.018 * scale, "distance": 0.10}
            else:
                res[name] = {"resolution": 0.045 * scale, "distance": 0.25}
    shapes["sub"] = box(0.0, y_lo, bx["xmax"], -0.30)
    shapes["air"] = box(0.0, y_lo, bx["xmax"], y_hi)
    res["sub"] = {"resolution": 0.10 * scale, "distance": 0.4}
    res["air"] = {"resolution": 0.10 * scale, "distance": 0.4}
    mesh = from_meshio(mesh_from_OrderedDict(shapes, res, default_resolution_max=0.12 * scale, filename=msh))
    return mesh, cells


def eps_region(scheme, state, side, xc, yc):
    """Проницаемость в точке (xc, yc) сечения стороны 1 (подвод) или 2 (фазовращатель)."""
    g = SCHEMES[scheme]
    eps_pcm = complex(N_PCM[state], 1e-5) ** 2
    if yc < -T_PCM:
        return N_SIO2 ** 2
    if side == 1:
        if g["shared"]:
            if yc < 0.0:
                return eps_pcm if xc < W_PCM / 2 else 1.0
            if yc < g["t_feed"] and xc < g["w_feed"] / 2:
                return N_SI ** 2
            return 1.0
        # подвод на подложке: кремний от -T_PCM до -T_PCM + t_feed, PCM под ним нет
        if yc < -T_PCM + g["t_feed"] and xc < g["w_feed"] / 2:
            return N_SI ** 2
        return 1.0
    # сторона 2: буфер, плёнка, гребень
    if yc < 0.0:
        return eps_pcm if xc < W_PCM / 2 else 1.0
    if xc < W_CORE / 2:
        if yc < D_AU:
            return EPS_AU
        if yc < D_AU + T_RIDGE:
            return N_SI ** 2
    return 1.0


def assign(basis0, cells, scheme, state, side):
    eps = basis0.zeros(dtype=complex)
    eps[basis0.get_dofs(elements="air")] = 1.0
    eps[basis0.get_dofs(elements="sub")] = N_SIO2 ** 2
    for name, xc, yc in cells:
        eps[basis0.get_dofs(elements=name)] = eps_region(scheme, state, side, xc, yc)
    return eps


class Side:
    """Набор мод одного сечения на квадратурных точках общей сетки."""

    def __init__(self, modes, dx, sigma):
        self.n = np.array([complex(m.n_eff) for m in modes])
        self.lam = self.n ** 2 * K0 ** 2
        order = np.argsort(np.abs(self.lam - sigma))
        self.n = self.n[order]
        self.modes = [modes[i] for i in order]
        self.dx = dx
        ex, ey, hx, hy = [], [], [], []
        self.tm = []
        for m in self.modes:
            (Ex, Ey), Ez = m.basis.interpolate(m.E)
            (Hx, Hy), Hz = m.basis.interpolate(m.H)
            ex.append(np.asarray(Ex).ravel().astype(np.complex64))
            ey.append(np.asarray(Ey).ravel().astype(np.complex64))
            hx.append(np.asarray(Hx).ravel().astype(np.complex64))
            hy.append(np.asarray(Hy).ravel().astype(np.complex64))
            e2 = np.abs(Ex) ** 2 + np.abs(Ey) ** 2 + np.abs(Ez) ** 2
            self.tm.append(float(np.sum(np.abs(Ey) ** 2 * dx) / np.sum(e2 * dx)))
        self.ex, self.ey, self.hx, self.hy = (np.array(v) for v in (ex, ey, hx, hy))
        self.tm = np.array(self.tm)
        w = np.asarray(dx).ravel().astype(np.complex64)
        self.w = w
        # Ветвь корня. femwell берёт главный корень beta = sqrt(lambda): для распространяющихся мод это
        # волна вперёд (Re beta > 0), но у затухающих мод с Im lambda < 0 главный корень имеет Im beta < 0,
        # то есть растёт вдоль +z - это волна, бегущая К стыку, а не от него. Для прямого базиса
        # нужна Im beta >= 0: у таких мод меняется знак beta, а с ним и знак H (E_t не меняется).
        self.flipped = np.zeros(len(self.n), dtype=bool)
        for k in range(len(self.n)):
            if self.n[k].imag < -1e-4:
                self.n[k] = -self.n[k]
                self.hx[k] *= -1
                self.hy[k] *= -1
                self.flipped[k] = True
        # Нормировка. femwell делит поля на корень из сопряжённой мощности Re int (e x h^*).z, которая
        # у затухающих мод равна нулю, и поля раздуваются; несопряжённая норма int (e x h).z у
        # комплексных мод (пары lambda, lambda^*) тоже нулевая. Безопасна норма int |e_t|^2 dS = 1;
        # формулы сшивания используют полные матрицы U и не требуют diag(U) = 1.
        l2 = np.sqrt(np.real((np.abs(self.ex) ** 2 + np.abs(self.ey) ** 2) @ w)).astype(np.complex64)
        for arr in (self.ex, self.ey, self.hx, self.hy):
            arr /= l2[:, None]
        # U[m,k] = int (e_m x h_k).z = int (ex_m hy_k - ey_m hx_k) - несопряжённая матрица
        self.U = (self.ex * w) @ self.hy.T - (self.ey * w) @ self.hx.T
        # C[m,k] = int (e_m x h_k^*).z  - сопряжённая (мощностная) матрица; Re C[m,m] - поток мощности моды
        self.C = (self.ex * w) @ np.conj(self.hy).T - (self.ey * w) @ np.conj(self.hx).T
        d = np.abs(np.diag(self.U))
        reg = d > 1e-3 * d.max()
        Un = self.U[np.ix_(reg, reg)] / np.sqrt(np.outer(d[reg], d[reg]))
        self.biorth = float(np.abs(Un - np.diag(np.diag(Un))).max())
        self.n_complex = int((~reg).sum())

    def cross(self, other):
        """X[n,m] = int (e_n(self) x h_m(other)).z"""
        return (self.ex * self.w) @ other.hy.T - (self.ey * self.w) @ other.hx.T


def solve_interface(Ua, Ub, X, Y, i, N, variant=1):
    """Падение моды i стороны a; усечение до первых N мод каждой стороны."""
    Ua, Ub = Ua[:N, :N], Ub[:N, :N]
    X, Y = X[:N, :N], Y[:N, :N]
    d = np.zeros(N, dtype=complex)
    d[i] = 1.0
    if variant == 1:
        UaTinv_XT = np.linalg.solve(Ua.T, X.T)          # U_a^{-T} X^T
        A = Ub + X @ UaTinv_XT
        t = np.linalg.solve(A, 2 * X @ d)
        r = UaTinv_XT @ t - d
    else:
        Uainv_Y = np.linalg.solve(Ua, Y)                # U_a^{-1} Y
        A = Ub.T + Y.T @ Uainv_Y
        t = np.linalg.solve(A, 2 * Y.T @ d)
        r = d - Uainv_Y @ t
    return t, r


def powers(t, r, Ca, Cb, N, ia):
    """Полные потоки прошедшего и отражённого полей, отнесённые к потоку падающей моды ia."""
    p_inc = Ca[ia, ia].real
    Ca, Cb = Ca[:N, :N], Cb[:N, :N]
    p_t = float(np.real(t @ Cb @ np.conj(t))) / p_inc
    p_r = float(np.real(r @ Ca @ np.conj(r))) / p_inc
    return p_t, p_r


def pick_feed(side):
    """Основная квази-TM мода подвода: без потерь, выше подложки, наибольшая доля |Ey|^2 среди верхних."""
    cand = [k for k in range(len(side.n)) if side.n[k].real > N_SIO2 and abs(side.n[k].imag) < 1e-3 and side.tm[k] > 0.5]
    if not cand:
        raise RuntimeError("нет направляемой квази-TM моды подвода")
    return max(cand, key=lambda k: side.n[k].real)


def pick_shifter(side, ref):
    """Плазмонная мода: ближайшая по комплексному показателю к опорному значению (продолжение корня)."""
    return int(np.argmin(np.abs(np.array([complex(z.real, abs(z.imag)) for z in side.n]) - ref)))


def main():
    scheme = sys.argv[1] if len(sys.argv) > 1 else "sync405"
    state = sys.argv[2] if len(sys.argv) > 2 else "a"
    boxname = sys.argv[3] if len(sys.argv) > 3 else "A"
    nmax = int(sys.argv[4]) if len(sys.argv) > 4 else 60
    scale = float(sys.argv[5]) if len(sys.argv) > 5 else 1.0
    n_guess = float(sys.argv[6]) if len(sys.argv) > 6 else 1.6
    tag = "%s_%s_%s_N%d_s%.2f_g%.2f" % (scheme, state, boxname, nmax, scale, n_guess)
    logf = os.path.join(HERE, "mm_status.txt")
    log("START %s %s" % (tag, time.strftime("%Y-%m-%d %H:%M:%S")), logf)
    from skfem import Basis, ElementTriP0
    from femwell.maxwell.waveguide import compute_modes

    t0 = time.time()
    msh = os.path.join(HERE, "mm_mesh_%s.msh" % tag)
    mesh, cells = build_mesh(scheme, boxname, scale, msh)
    basis0 = Basis(mesh, ElementTriP0())
    log("  сетка: %d элементов, %.0f с" % (mesh.t.shape[1], time.time() - t0), logf)

    sides = {}
    sigma = K0 ** 2 * n_guess ** 2
    for sd in (1, 2):
        eps = assign(basis0, cells, scheme, state, sd)
        t1 = time.time()
        modes = compute_modes(basis0, eps, wavelength=LAM, num_modes=nmax, order=1, n_guess=n_guess)
        dx = modes[0].basis.dx
        sides[sd] = Side(list(modes), dx, sigma)
        log("  сторона %d: %d мод за %.0f с; квадратурных точек %d" % (sd, nmax, time.time() - t1, sides[sd].ex.shape[1]), logf)
        # проверка биортогональности внутри стороны: наибольший внедиагональный элемент U при diag(U) = 1
        log("  сторона %d: биортогональность обычных мод max|U_mk|/sqrt(U_mm U_kk) = %.2e; комплексных мод (U_mm ~ 0): %d; "
            "ветвь переключена у %d; мод с Re n > %.3f: %d, затухающих (Re n < 0.05): %d"
            % (sd, sides[sd].biorth, sides[sd].n_complex, int(sides[sd].flipped.sum()), N_SIO2,
               int(np.sum(sides[sd].n.real > N_SIO2)), int(np.sum(np.abs(sides[sd].n.real) < 0.05))), logf)
    s1, s2 = sides[1], sides[2]
    g1 = pick_feed(s1)
    g2 = pick_shifter(s2, REF_SHIFTER[state])
    log("  подвод: n = %.6f%+.2ei (TM %.3f, ранг %d); фазовращатель: n = %.6f%+.5fi (TM %.3f, ранг %d)"
        % (s1.n[g1].real, s1.n[g1].imag, s1.tm[g1], g1, s2.n[g2].real, s2.n[g2].imag, s2.tm[g2], g2), logf)

    # спектры для рисунка
    with io.open(os.path.join(HERE, "mm_spectrum_%s.csv" % tag), "w", encoding="utf-8", newline="\n") as f:
        f.write("side,rank,re_neff,im_neff,tm_fraction,guided\n")
        for sd, s, gi in ((1, s1, g1), (2, s2, g2)):
            for k in range(len(s.n)):
                f.write("%d,%d,%.7f,%.4e,%.4f,%d\n" % (sd, k, s.n[k].real, s.n[k].imag, s.tm[k], 1 if k == gi else 0))

    X12 = s2.cross(s1)      # X[n,m] = int e_2n x h_1m  (сторона a = 1, b = 2)
    Y12 = s1.cross(s2)      # Y[m,n] = int e_1m x h_2n
    # первый порядок - перекрытие по потоку (для сравнения с прежними расчётами)
    c12 = Y12[g1, g2]
    c21 = X12[g2, g1]
    Cc = s1.ex * s1.w @ np.conj(s2.hy).T - s1.ey * s1.w @ np.conj(s2.hx).T   # int e_1 x h_2^*
    Cc2 = s2.ex * s2.w @ np.conj(s1.hy).T - s2.ey * s2.w @ np.conj(s1.hx).T  # int e_2 x h_1^*
    eta_conj = abs(Cc[g1, g2] + np.conj(Cc2[g2, g1])) ** 2 / (4 * abs(s1.C[g1, g1]) * abs(s2.C[g2, g2]))
    eta_unconj = abs(c12 * c21 / (s1.U[g1, g1] * s2.U[g2, g2]))
    log("  первый порядок: eta_conj = %.5f (%.3f дБ), eta_unconj = %.5f" % (eta_conj, -10 * math.log10(eta_conj), eta_unconj), logf)
    # rho = |N|/P - отношение несопряжённой нормы моды к её потоку мощности: единица у моды без потерь,
    # меньше единицы у плазмонной; по взаимности T_rev/T_fwd = (rho2/rho1)^2
    rho1 = abs(s1.U[g1, g1]) / s1.C[g1, g1].real
    rho2 = abs(s2.U[g2, g2]) / s2.C[g2, g2].real
    log("  rho = |N|/P: подвод %.6f, фазовращатель %.6f; предсказание взаимности T_rev/T_fwd = (rho2/rho1)^2 = %.5f"
        % (rho1, rho2, (rho2 / rho1) ** 2), logf)

    out = io.open(os.path.join(HERE, "mm_results.csv"), "a", encoding="utf-8", newline="\n")
    if out.tell() == 0:
        out.write("tag,scheme,state,box,scale,elements,N,dir,variant,T_guided,R_guided,P_trans,P_refl,balance,"
                  "t_re,t_im,r_re,r_im,eta_conj,eta_unconj,n1_re,n1_im,n2_re,n2_im,rho1,rho2,biorth1,biorth2\n")
    Ns = sorted(set([n for n in (5, 10, 15, 20, 30, 40, 50, 60, 80, 100, 120, 150, 200, 250, 300) if n <= nmax] + [nmax]))
    last = {}
    for N in Ns:
        for direction, (sa, sb, X, Y, ia, ib) in (("12", (s1, s2, X12, Y12, g1, g2)), ("21", (s2, s1, Y12.T, X12.T, g2, g1))):
            # для направления 2->1: сторона a = 2, b = 1; X[n,m] = int e_1n x h_2m = Y12[n,m]... см. док-строку:
            # X_{ba}[n,m] = int e_bn x h_am. При a=2,b=1: int e_1n x h_2m = Y12[n,m], то есть X = Y12; Y = X12^T? нет:
            # Y_{ab}[m,n] = int e_am x h_bn = int e_2m x h_1n = X12[m,n]. Поэтому для 21: X = Y12, Y = X12.
            if direction == "21":
                X, Y = Y12, X12
            if ia >= N or ib >= N:
                continue
            for variant in (1, 2):
                t, r = solve_interface(sa.U, sb.U, X, Y, ia, N, variant)
                p_t, p_r = powers(t, r, sa.C, sb.C, N, ia)
                T = abs(t[ib]) ** 2 * sb.C[ib, ib].real / sa.C[ia, ia].real
                R = abs(r[ia]) ** 2
                bal = 1.0 - p_t - p_r
                out.write("%s,%s,%s,%s,%.2f,%d,%d,%s,%d,%.6f,%.6f,%.6f,%.6f,%.2e,%.6e,%.6e,%.6e,%.6e,%.6f,%.6f,%.7f,%.4e,%.7f,%.5e,%.6f,%.3e,%.6f,%.3e\n"
                          % (tag, scheme, state, boxname, scale, mesh.t.shape[1], N, direction, variant, T, R, p_t, p_r, bal,
                             t[ib].real, t[ib].imag, r[ia].real, r[ia].imag, eta_conj, eta_unconj,
                             s1.n[g1].real, s1.n[g1].imag, s2.n[g2].real, s2.n[g2].imag,
                             rho1, rho2, s1.biorth, s2.biorth))
                last[(direction, variant)] = (t, r, T, R, p_t, p_r, bal)
                if N == nmax or variant == 1:
                    log("  N=%3d %s v%d: T=%.5f (%.3f дБ) R=%.2e P_t=%.5f P_r=%.5f баланс %.1e"
                        % (N, direction, variant, T, -10 * math.log10(T), R, p_t, p_r, bal), logf)
    out.close()

    # взаимность: t_12 N_2 = t_21 N_1 при несопряжённых нормах N = U_gg (проверка на полном базисе)
    N1, N2 = s1.U[g1, g1], s2.U[g2, g2]
    for variant in (1, 2):
        t12 = last[("12", variant)][0][g2]
        t21 = last[("21", variant)][0][g1]
        log("  взаимность (вариант %d): t12*N2 = %.6f%+.6fi, t21*N1 = %.6f%+.6fi, |отн. разность| = %.2e"
            % (variant, (t12 * N2).real, (t12 * N2).imag, (t21 * N1).real, (t21 * N1).imag,
               abs(t12 * N2 - t21 * N1) / abs(t12 * N2)), logf)
    log("  |N1/N2|^2 = %.5f: во столько раз передача по мощности плазмон -> диэлектрик больше обратной"
        % (abs(N1 / N2) ** 2 * (s1.C[g1, g1].real / s2.C[g2, g2].real) ** 0), logf)

    # поля на плоскости стыка для рисунка: падающая мода, прошедшее поле (сторона 2), отражённое (сторона 1)
    t, r = last[("12", 1)][0], last[("12", 1)][1]
    N = nmax
    xy = s1.modes[0].basis.global_coordinates()
    x = np.asarray(xy[0]).ravel()
    y = np.asarray(xy[1]).ravel()
    ey_inc = s1.ey[g1]
    ey_tr = (t[:N, None] * s2.ey[:N]).sum(axis=0)
    ey_ref = (r[:N, None] * s1.ey[:N]).sum(axis=0)
    ey_g2 = s2.ey[g2]
    np.savez_compressed(os.path.join(HERE, "mm_fields_%s.npz" % tag), x=x, y=y, ey_inc=ey_inc, ey_tr=ey_tr,
                        ey_ref=ey_ref, ey_g2=ey_g2, t=t, r=r, n1=s1.n, n2=s2.n,
                        eps1=np.asarray(basis0.interpolate(assign(basis0, cells, scheme, state, 1))).ravel(),
                        eps2=np.asarray(basis0.interpolate(assign(basis0, cells, scheme, state, 2))).ravel())
    try:
        os.remove(msh)
    except OSError:
        pass
    log("DONE %s за %.0f с" % (tag, time.time() - t0), logf)


if __name__ == "__main__":
    main()
