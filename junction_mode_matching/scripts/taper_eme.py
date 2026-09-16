"""Переход по ширине золота как каскад стыков: разложение по местным модам (EME) для любой формы острия.

Базис - taper_basis.py: местные моды сечений на сетке ширин золота 0...450 нм плюс подводящий гребень,
все на одной сетке, с матрицами Грама U (несопряжённая), C (сопряжённая) и перекрёстными X, Y между
соседними сечениями. Профиль острия w(z) монотонен, поэтому лестница из ширин сетки проходит все стыки
по порядку: S-матрица перехода = каскад S-матриц стыков и участков распространения.

S-матрица стыка a|b (нормировка мод внутри каждого сечения общая для обоих его стыков):
    T_ba = (U_b + X U_a^{-T} X^T)^{-1} 2X,   R_aa = U_a^{-T} X^T T_ba - I    (падение слева, для всех мод)
    T_ab = (U_a + Y U_b^{-T} Y^T)^{-1} 2Y,   R_bb = U_b^{-T} Y^T T_ab - I    (падение справа)
где X[n,m] = int e_bn x h_am, Y[m,n] = int e_am x h_bn. Участок длиной d: P = diag(exp(i beta d)),
Im beta >= 0 у всех мод базиса, так что затухающие моды убывают от обоих стыков. Каскад - формула
Редхеффера, устойчивая при затухающих модах.

Мера: доля мощности, дошедшая в направляемую моду фазовращателя (вперёд) или подводящего гребня (назад),
отнесённая к потоку падающей моды; поглощение вдоль перехода сидит в комплексных beta местных мод.
Опорная кривая - идеально адиабатический проход exp(-int alpha(w(z)) dz) по прослеженной ветви.

Запуск: python taper_eme.py <basis.npz> [N] [профили через запятую] [длины через запятую, мкм]
"""
import io
import math
import os
import sys

import numpy as np

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

HERE = os.path.dirname(os.path.abspath(__file__))
K0 = 2 * math.pi / 1.55
W_FULL = 450.0


def profile(shape, u, table=None):
    """Полуширина металла в долях рабочей как функция u = z/L (те же формулы, что в RunTaper.java)."""
    u = np.clip(u, 0.0, 1.0)
    if shape == "linear":
        return u
    if shape == "sqrt":
        return np.sqrt(u)
    if shape == "quad":
        return u ** 2
    if shape == "quart":
        return u ** 0.25
    if shape == "cubic":
        return u ** 3
    if shape == "ellipse":
        return np.sqrt(1.0 - (1.0 - u) ** 2)
    if shape == "ellipsec":
        return 1.0 - np.sqrt(1.0 - u ** 2)
    if shape == "rcos":
        return 0.5 * (1.0 - np.cos(np.pi * u))
    if shape == "expo":
        a = 3.0
        return (np.exp(a * u) - 1.0) / (np.exp(a) - 1.0)
    if shape == "gauss":
        s = 0.35
        g0 = math.exp(-1.0 / (2 * s * s))
        return (np.exp(-(1.0 - u) ** 2 / (2 * s * s)) - g0) / (1.0 - g0)
    if shape == "klop":
        _init_klop()
        return np.interp(u, _KLOP_GRID[0], _KLOP_GRID[1])
    if shape == "fastwin":
        return _piecewise(u, 0.45, 0.50)
    if shape == "slowwin":
        return _piecewise(u, 0.25, 0.75)
    if shape == "table":
        t = np.asarray(table)
        return np.interp(u, t[:, 0], t[:, 1])
    raise ValueError(shape)


def _piecewise(u, u0, u1):
    f_lo, f_hi = 210.0 / 450.0, 290.0 / 450.0
    out = np.where(u <= u0, f_lo * u / u0, np.where(u <= u1, f_lo + (f_hi - f_lo) * (u - u0) / (u1 - u0),
                                                       f_hi + (1.0 - f_hi) * (u - u1) / (1.0 - u1)))
    return out


def _klop(u):
    a = 4.0

    def phi(x):
        n = 200
        s = 0.0
        for i in range(n):
            y0, y1 = x * i / n, x * (i + 1) / n
            s += 0.5 * (kern(y0) + kern(y1)) * (y1 - y0)
        return s

    def kern(y):
        r = 1.0 - y * y
        if r <= 0:
            return 0.5
        t = a * math.sqrt(r)
        return _i1(t) / t

    raw = 0.5 + phi(2 * u - 1) / phi(1.0) * 0.5
    lo = 0.5 + phi(-1.0) / phi(1.0) * 0.5
    hi = 0.5 + phi(1.0) / phi(1.0) * 0.5
    return (raw - lo) / (hi - lo)


def _i1(x):
    term = x / 2.0
    s = term
    for k in range(1, 30):
        term *= (x * x / 4.0) / (k * (k + 1.0))
        s += term
    return s


_KLOP_GRID = None


def _init_klop():
    global _KLOP_GRID
    if _KLOP_GRID is None:
        ug = np.linspace(0.0, 1.0, 801)
        _KLOP_GRID = (ug, np.array([_klop(float(x)) for x in ug]))


_PROFILE_CACHE = {}


def z_of_width(shape, w_nm, L, table=None):
    """Координата z, на которой профиль достигает полуширины w (мкм), обратной функцией по таблице."""
    key = (shape, None if table is None else tuple(map(tuple, np.asarray(table))))
    if key not in _PROFILE_CACHE:
        if shape == "klop":
            _init_klop()
        u = np.linspace(0.0, 1.0, 20001)
        _PROFILE_CACHE[key] = (u, profile(shape, u, table) * W_FULL)
    u, h = _PROFILE_CACHE[key]
    # профиль монотонен неубывающий; берём первое достижение
    return float(np.interp(w_nm, h, u * L))


class Basis:
    def __init__(self, path, N):
        d = np.load(path, allow_pickle=True)
        self.labels = [str(x) for x in d["labels"]]
        # ширины из имён сечений, чтобы читать и частично записанный базис
        self.widths = [-1 if lab == "feed" else int(lab[1:]) for lab in self.labels]
        self.N = N
        self.n = {lab: d["n_" + lab][:N] for lab in self.labels}
        self.U = {lab: d["U_" + lab][:N, :N] for lab in self.labels}
        self.C = {lab: d["C_" + lab][:N, :N] for lab in self.labels}
        self.g = {lab: int(d["g_" + lab]) for lab in self.labels}
        self.X = {}
        self.Y = {}
        for a, b in zip(self.labels[:-1], self.labels[1:]):
            self.X[(a, b)] = d["X_%s_%s" % (a, b)][:N, :N]
            self.Y[(a, b)] = d["Y_%s_%s" % (a, b)][:N, :N]
        self._iface = {}

    def interface(self, a, b):
        """Четыре блока S-матрицы стыка a|b: T_ba, R_aa, T_ab, R_bb (считаются один раз на стык)."""
        if (a, b) in self._iface:
            return self._iface[(a, b)]
        self._iface[(a, b)] = self._interface(a, b)
        return self._iface[(a, b)]

    def _interface(self, a, b):
        Ua, Ub, X, Y = self.U[a], self.U[b], self.X[(a, b)], self.Y[(a, b)]
        N = self.N
        I = np.eye(N, dtype=complex)
        # падение слева: непрерывность E свёрнута с h_a, непрерывность H - с e_b
        UaTinv_XT = np.linalg.solve(Ua.T, X.T)          # U_a^{-T} X^T
        T_ba = np.linalg.solve(Ub + X @ UaTinv_XT, 2 * X)
        R_aa = UaTinv_XT @ T_ba - I
        # падение справа: те же формулы с переставленными сторонами, перекрёстная матрица - Y[m,n] = int e_am x h_bn
        UbTinv_YT = np.linalg.solve(Ub.T, Y.T)          # U_b^{-T} Y^T
        T_ab = np.linalg.solve(Ua + Y @ UbTinv_YT, 2 * Y)
        R_bb = UbTinv_YT @ T_ab - I
        return T_ba, R_aa, T_ab, R_bb


def redheffer(S1, S2):
    """Каскад двух S-матриц (T_ba, R_aa, T_ab, R_bb): сначала S1, затем S2."""
    T1, R1, T1b, R1b = S1
    T2, R2, T2b, R2b = S2
    N = T1.shape[0]
    I = np.eye(N, dtype=complex)
    M = np.linalg.solve(I - R1b @ R2, I)          # (I - R1b R2)^{-1}
    T_ba = T2 @ M @ T1
    R_aa = R1 + T1b @ R2 @ M @ T1
    Mb = np.linalg.solve(I - R2 @ R1b, I)
    T_ab = T1b @ Mb @ T2b
    R_bb = R2b + T2 @ R1b @ Mb @ T2b
    return T_ba, R_aa, T_ab, R_bb


def propagation(n, d):
    """Участок однородного волновода длиной d (мкм): моды бегут с exp(i beta d), Im beta >= 0."""
    beta = K0 * n
    P = np.diag(np.exp(1j * beta * d))
    Z = np.zeros_like(P)
    return P, Z, P, Z


def taper_smatrix(basis, shape, L, table=None, l_in=0.0, l_out=0.0):
    """S-матрица перехода: подводящий гребень | сердцевина без металла | лестница ширин | фазовращатель.

    Участок сетки с шириной w_i занимает z от середины отрезка (w_{i-1}, w_i) до середины (w_i, w_{i+1});
    первая ширина (0) начинается при z = 0, последняя (450) тянется до конца перехода и дальше на l_out.
    """
    labs = basis.labels
    widths = basis.widths
    # границы участков по z
    zs = []
    for i in range(1, len(labs)):          # labs[0] = feed
        w = widths[i]
        if i == 1:
            z0 = 0.0
        else:
            z0 = z_of_width(shape, 0.5 * (widths[i - 1] + w), L, table)
        if i == len(labs) - 1:
            z1 = L + l_out
        else:
            z1 = z_of_width(shape, 0.5 * (w + widths[i + 1]), L, table)
        zs.append((z0, max(z1, z0)))
    S = propagation(basis.n[labs[0]], l_in)
    S = redheffer(S, basis.interface(labs[0], labs[1]))
    for i in range(1, len(labs)):
        z0, z1 = zs[i - 1]
        S = redheffer(S, propagation(basis.n[labs[i]], z1 - z0))
        if i < len(labs) - 1:
            S = redheffer(S, basis.interface(labs[i], labs[i + 1]))
    return S, zs


def adiabatic_reference(basis, shape, L, table=None):
    """exp(-int alpha dz) по прослеженной ветви, alpha_i на участке ширины w_i."""
    labs, widths = basis.labels, basis.widths
    total = 0.0
    for i in range(1, len(labs)):
        w = widths[i]
        z0 = 0.0 if i == 1 else z_of_width(shape, 0.5 * (widths[i - 1] + w), L, table)
        z1 = L if i == len(labs) - 1 else z_of_width(shape, 0.5 * (w + widths[i + 1]), L, table)
        g = basis.g[labs[i]]
        alpha = 2 * K0 * abs(basis.n[labs[i]][g].imag)
        total += alpha * max(z1 - z0, 0.0)
    return math.exp(-total)


def evaluate(basis, shape, L, table=None):
    S, zs = taper_smatrix(basis, shape, L, table)
    T_ba, R_aa, T_ab, R_bb = S
    a, b = basis.labels[0], basis.labels[-1]
    ga, gb = basis.g[a], basis.g[b]
    Ca, Cb = basis.C[a], basis.C[b]
    # вперёд: подводящий гребень -> фазовращатель; амплитуды в конце перехода (z = L), без хвоста
    t_fwd = T_ba[:, ga]
    T_fwd = abs(t_fwd[gb]) ** 2 * Cb[gb, gb].real / Ca[ga, ga].real
    P_fwd = float(np.real(t_fwd @ Cb @ np.conj(t_fwd))) / Ca[ga, ga].real
    r_fwd = R_aa[:, ga]
    R_fwd = abs(r_fwd[ga]) ** 2
    # назад: из плазмонной моды в моду гребня
    t_rev = T_ab[:, gb]
    T_rev = abs(t_rev[ga]) ** 2 * Ca[ga, ga].real / Cb[gb, gb].real
    P_rev = float(np.real(t_rev @ Ca @ np.conj(t_rev))) / Cb[gb, gb].real
    r_rev = R_bb[:, gb]
    R_rev = abs(r_rev[gb]) ** 2
    T_ad = adiabatic_reference(basis, shape, L, table)
    # мера трёхмерной модели: |S21|^2 делится на затухание рабочей моды по ВСЕЙ длине от стыка до порта,
    # то есть участку перехода зачитывается затухание рабочей моды; та же мера здесь - множитель exp(alpha_450 L)
    alpha_end = 2 * K0 * abs(basis.n[b][gb].imag)
    eq = math.exp(alpha_end * L)
    return dict(T_fwd=T_fwd, T_rev=T_rev, P_fwd=P_fwd, P_rev=P_rev, R_fwd=R_fwd, R_rev=R_rev, T_ad=T_ad,
                T_fwd_eq=T_fwd * eq, T_rev_eq=T_rev * eq, T_ad_eq=T_ad * eq)


def main():
    path = sys.argv[1]
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 100
    shapes = sys.argv[3].split(",") if len(sys.argv) > 3 else ["step", "linear", "quad", "sqrt", "ellipse", "ellipsec", "expo", "gauss", "rcos", "klop", "fastwin", "slowwin"]
    lengths = [float(x) for x in sys.argv[4].split(",")] if len(sys.argv) > 4 else [0.25, 0.5, 1.0, 2.0]
    basis = Basis(path, N)
    direct = load_direct(basis, path.replace("taper_basis_", "taper_basis_direct_"))
    print("базис: %d сечений, N = %d, ширины %s; прямые стыки: %s" % (len(basis.labels), N, basis.widths, "есть" if direct is not None else "нет"))
    print("ветвь: " + " ".join("%d:%.4f%+.4fi" % (w, basis.n[l][basis.g[l]].real, abs(basis.n[l][basis.g[l]].imag)) for w, l in zip(basis.widths, basis.labels)))
    out = io.open(os.path.join(HERE, "taper_eme_results.csv"), "a", encoding="utf-8", newline="\n")
    if out.tell() == 0:
        out.write("basis,N,shape,L_um,T_fwd,T_rev,P_fwd,P_rev,R_fwd,R_rev,T_ad,loss_fwd_dB,loss_rev_dB,loss_ad_dB,conv_fwd_dB,T_fwd_eq,T_rev_eq,T_ad_eq\n")
    for L in lengths:
        for sh in shapes:
            if sh == "step":
                # контроль: металл начинается скачком при z = L, то есть стык 0 -> 450 после участка без металла
                res = evaluate_step(basis, L, direct)
            elif sh == "butt":
                # контроль: соединение встык sync405 без участка без металла (прямой стык гребень -> плёнка)
                if direct is None:
                    continue
                res = evaluate_butt(basis, direct)
            else:
                res = evaluate(basis, sh, L)
            lf, lr, la = (-10 * math.log10(res["T_fwd"]), -10 * math.log10(res["T_rev"]), -10 * math.log10(res["T_ad"]))
            print("L=%.2f %-9s T_fwd=%.4f (%.3f дБ) T_rev=%.4f (%.3f дБ) адиаб. %.4f (%.3f дБ) преобр. %.3f дБ  R_fwd=%.1e R_rev=%.1e  P_fwd=%.4f P_rev=%.4f | в мере COMSOL %.4f / %.4f, адиаб. %.4f"
                  % (L, sh, res["T_fwd"], lf, res["T_rev"], lr, res["T_ad"], la, lf - la, res["R_fwd"], res["R_rev"], res["P_fwd"], res["P_rev"],
                     res["T_fwd_eq"], res["T_rev_eq"], res["T_ad_eq"]))
            out.write("%s,%d,%s,%.3f,%.6f,%.6f,%.6f,%.6f,%.3e,%.3e,%.6f,%.4f,%.4f,%.4f,%.4f,%.6f,%.6f,%.6f\n"
                      % (os.path.basename(path), N, sh, L, res["T_fwd"], res["T_rev"], res["P_fwd"], res["P_rev"], res["R_fwd"], res["R_rev"], res["T_ad"], lf, lr, la, lf - la,
                         res["T_fwd_eq"], res["T_rev_eq"], res["T_ad_eq"]))
            out.flush()
    out.close()


def load_direct(basis, path):
    """Самостоятельный набор для контролей: моды feed, w000, w450 одного прогона и их стыки (taper_basis_direct.py).

    Возвращает объект с полями labels, n, U, C, g, X, Y и методом interface - тот же интерфейс, что у Basis,
    поэтому evaluate_butt и evaluate_step работают только внутри этого набора и не смешивают его с основным.
    """
    if not os.path.exists(path):
        return None
    d = np.load(path, allow_pickle=True)
    if "U_feed" not in d:
        return None
    N = basis.N

    class Direct:
        pass

    o = Direct()
    o.labels = ["feed", "w000", "w450"]
    o.N = N
    o.n = {lab: d["n_" + lab][:N] for lab in o.labels}
    o.U = {lab: d["U_" + lab][:N, :N] for lab in o.labels}
    o.C = {lab: d["C_" + lab][:N, :N] for lab in o.labels}
    o.g = {lab: int(d["g_" + lab]) for lab in o.labels}
    o.X = {}
    o.Y = {}
    for a, b in (("feed", "w000"), ("w000", "w450"), ("feed", "w450")):
        o.X[(a, b)] = d["X_%s_%s" % (a, b)][:N, :N]
        o.Y[(a, b)] = d["Y_%s_%s" % (a, b)][:N, :N]
    o._iface = {}
    o.interface = lambda a, b: Basis.interface(o, a, b)
    o._interface = lambda a, b: Basis._interface(o, a, b)
    return o


def evaluate_butt(basis, direct):
    """Соединение встык sync405: подводящий гребень сразу в полную плёнку, без участка без металла."""
    a, b = "feed", "w450"
    T_ba, R_aa, T_ab, R_bb = direct.interface(a, b)
    ga, gb = direct.g[a], direct.g[b]
    Ca, Cb = direct.C[a], direct.C[b]
    t_fwd, t_rev = T_ba[:, ga], T_ab[:, gb]
    T_fwd = abs(t_fwd[gb]) ** 2 * Cb[gb, gb].real / Ca[ga, ga].real
    T_rev = abs(t_rev[ga]) ** 2 * Ca[ga, ga].real / Cb[gb, gb].real
    return dict(T_fwd=T_fwd, T_rev=T_rev,
                P_fwd=float(np.real(t_fwd @ Cb @ np.conj(t_fwd))) / Ca[ga, ga].real,
                P_rev=float(np.real(t_rev @ Ca @ np.conj(t_rev))) / Cb[gb, gb].real,
                R_fwd=abs(R_aa[ga, ga]) ** 2, R_rev=abs(R_bb[gb, gb]) ** 2, T_ad=1.0,
                T_fwd_eq=T_fwd, T_rev_eq=T_rev, T_ad_eq=1.0)


def evaluate_step(basis, L, direct=None):
    """Контроль без перехода: сердцевина без металла длиной L, затем скачок к полной ширине.

    Прямой стык 0 -> 450 берётся из taper_basis_direct.npz (перекрёстные матрицы той пары на той же
    сетке); без него стык собирается лестницей стыков нулевой длины, что в усечённом базисе неточно.
    """
    if direct is not None:
        # целиком на наборе прямого прогона: гребень | сердцевина без металла длиной L | скачок к плёнке
        src = direct
        S = direct.interface("feed", "w000")
        S = redheffer(S, propagation(direct.n["w000"], L))
        S = redheffer(S, direct.interface("w000", "w450"))
        a, b = "feed", "w450"
    else:
        src = basis
        labs = basis.labels
        S = basis.interface(labs[0], labs[1])
        S = redheffer(S, propagation(basis.n[labs[1]], L))
        for i in range(1, len(labs) - 1):
            S = redheffer(S, basis.interface(labs[i], labs[i + 1]))
        a, b = labs[0], labs[-1]
    T_ba, R_aa, T_ab, R_bb = S
    ga, gb = src.g[a], src.g[b]
    Ca, Cb = src.C[a], src.C[b]
    t_fwd = T_ba[:, ga]
    t_rev = T_ab[:, gb]
    alpha_end = 2 * K0 * abs(src.n[b][gb].imag)
    eq = math.exp(alpha_end * L)
    T_fwd = abs(t_fwd[gb]) ** 2 * Cb[gb, gb].real / Ca[ga, ga].real
    T_rev = abs(t_rev[ga]) ** 2 * Ca[ga, ga].real / Cb[gb, gb].real
    return dict(T_fwd=T_fwd, T_rev=T_rev,
                P_fwd=float(np.real(t_fwd @ Cb @ np.conj(t_fwd))) / Ca[ga, ga].real,
                P_rev=float(np.real(t_rev @ Ca @ np.conj(t_rev))) / Cb[gb, gb].real,
                R_fwd=abs(R_aa[ga, ga]) ** 2, R_rev=abs(R_bb[gb, gb]) ** 2, T_ad=1.0,
                T_fwd_eq=T_fwd * eq, T_rev_eq=T_rev * eq, T_ad_eq=eq)


if __name__ == "__main__":
    main()
