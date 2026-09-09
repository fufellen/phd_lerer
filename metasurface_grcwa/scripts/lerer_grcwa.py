"""Общий модуль расчётов метаповерхностей А.М. Лерера открытым решателем grcwa (метод фурье-мод).

Что здесь есть:
- материалы в оптической конвенции exp(-i w t), Im eps >= 0 для пассивной среды:
  золото по Джонсону-Кристи (таблица composite_ema, 187.9-1937 нм), золото по таблице Лерера
  Au.txt (200-2000 нм, шаг 1 нм) и продолжение обеих таблиц в ИК моделью Друде с плазменной
  частотой 9.03 эВ и затуханием, подобранным по непрерывности на стыке;
- правила смешивания композита «наночастицы в матрице»: полный Максвелл-Гарнетт (у Лерера -
  «формула Клаузиуса-Моссотти»), первый порядок по концентрации, формула Лерера, Бруггеман,
  Лоренц-Лоренц, Максвелл-Гарнетт с поправкой MLWA к поляризуемости сферы, Максвелл-Гарнетт для
  нескольких сортов частиц (исправленная COMPOSITE_3);
- сборка двумерно-периодической стопки слоёв в grcwa: однородные слои и слои с столбиками
  квадратного, круглого и эллиптического сечения; конус задаётся лестницей из тонких слоёв;
- поглощение по слоям через объёмный интеграл Im eps |E|^2 решателя и среднее |E|^2 в области
  композита - то, что нужно самосогласованному нелинейному расчёту.

Конвенции. Длины в grcwa - микрометры, freq = 1/lambda, скорость света 1. Таблицы composite_ema
хранят eps в конвенции Лерера (Im < 0) и здесь комплексно сопрягаются. Решатель grcwa - в
оптической конвенции (Im eps > 0), как в rcwa_external_check.py главы 4 диссертации и в
absorber_rcwa.py заметки об улучшении поглотителя, где он сверен с COMSOL (расхождение 0.007).
"""
from __future__ import annotations

import cmath
import io
import math
import os
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO / "composite_ema" / "scripts"))
import composite_ema as ce  # noqa: E402

import grcwa  # noqa: E402
import grcwa.rcwa as _rc  # noqa: E402
from grcwa.fft_funs import get_conv  # noqa: E402
from scipy.optimize import minimize  # noqa: E402


def _matrix_zintegral_robust(q, thickness, shift=1e-12):
    """Замена grcwa.rcwa.Matrix_zintegral, устойчивая к вырожденным собственным числам.

    В однородном слое квадратной решётки порядки (1,0), (-1,0), (0,1), (0,-1) имеют одно и то же q,
    и исходная формула делит на q_j - conj(q_i) = 0 вне диагонали (сдвиг shift стоит только на
    диагонали) - объёмный интеграл возвращает NaN. Здесь для |q_j - conj(q_i)| -> 0 и
    |q_j + conj(q_i)| -> 0 подставлены пределы t и t exp(i q_j t).
    """
    qi, qj = _rc.Gmeshgrid(q)
    qij = qj - np.conj(qi)
    small = np.abs(qij) < 1e-9
    qij_safe = np.where(small, 1.0, qij)
    Maa = np.where(small, thickness + 0j, (np.exp(1j * qij_safe * thickness) - 1.0) / 1j / qij_safe)
    qij2 = qj + np.conj(qi)
    small2 = np.abs(qij2) < 1e-9
    qij2_safe = np.where(small2, 1.0, qij2)
    Mab = np.where(small2, thickness * np.exp(1j * qj * thickness),
                   (np.exp(1j * qj * thickness) - np.exp(-1j * np.conj(qi) * thickness)) / 1j / qij2_safe)
    tmp1 = np.vstack((Maa, Mab))
    tmp2 = np.vstack((Mab, Maa))
    return np.hstack((tmp1, tmp2))


_rc.Matrix_zintegral = _matrix_zintegral_robust

RESULTS = HERE.parent / "results"
DATA = HERE.parent / "data"
RESULTS.mkdir(parents=True, exist_ok=True)

EPS0 = 8.8541878128e-12
C_LIGHT = 2.99792458e8
HBAR_EV_S = 6.582119569e-16
HW_P_AU_EV = 9.03          # плазменная частота золота (Rakic 1998, Ordal 1985)
HGAMMA_ORDAL_EV = 0.0267   # затухание Друде по Ordal 1985 (215 см^-1)

N_HOST = 1.77
N_SUB = 1.45
C_FILL = 0.10

# ---------------------------------------------------------------------------
# Материалы (оптическая конвенция exp(-i w t), Im eps >= 0)
# ---------------------------------------------------------------------------

def eps_au_jc(lam_nm: float) -> complex:
    """Au по Johnson-Christy 1972 из таблицы composite_ema, сопряжено в оптическую конвенцию."""
    return ce.eps_au(float(lam_nm)).conjugate()


def eps_cu_jc(lam_nm: float) -> complex:
    return ce.eps_cu(float(lam_nm)).conjugate()


def eps_ag_jc(lam_nm: float) -> complex:
    return ce.eps_ag(float(lam_nm)).conjugate()


def _load_lerer_nk(path: Path = DATA / "Au_lerer_nk.txt"):
    lam, n, k = [], [], []
    for line in io.open(path, encoding="utf-8", errors="ignore"):
        p = line.split()
        if len(p) < 3:
            continue
        try:
            lam.append(float(p[0])); n.append(float(p[1])); k.append(float(p[2]))
        except ValueError:
            pass
    return np.array(lam), np.array(n), np.array(k)


_LER_LAM, _LER_N, _LER_K = _load_lerer_nk()
LERER_RANGE = (float(_LER_LAM[0]), float(_LER_LAM[-1]))
JC_RANGE = ce.MATERIAL_RANGE_NM[4]


def eps_au_lerer(lam_nm: float) -> complex:
    """Au по таблице Лерера Au.txt (n, k), линейная интерполяция по длине волны."""
    n = float(np.interp(lam_nm, _LER_LAM, _LER_N))
    k = float(np.interp(lam_nm, _LER_LAM, _LER_K))
    return complex(n, k) ** 2


def _drude(lam_nm: float, eps_inf: float, hw_p: float, hgamma: float) -> complex:
    hw = 1239.841984 / lam_nm
    return eps_inf - hw_p ** 2 / (hw ** 2 + 1j * hgamma * hw)


def _drude_fit(fn, lam_lo: float, lam_hi: float, eps_inf: float = 1.0) -> tuple[float, float, float]:
    """Плазменная частота и затухание Друде (эВ) по наименьшим квадратам к таблице на [lam_lo, lam_hi];
    eps_inf зафиксировано (1 - как у Ordal 1985). Возвращает (hw_p, hgamma, относительный скачок на стыке)."""
    lams = np.linspace(lam_lo, lam_hi, 25)
    tab = np.array([fn(l) for l in lams])

    def cost(p):
        d = np.array([_drude(l, eps_inf, p[0], p[1]) for l in lams])
        return float(np.sum(np.abs((d - tab) / tab) ** 2))

    res = minimize(cost, x0=[HW_P_AU_EV, 0.05], method="Nelder-Mead", options=dict(xatol=1e-6, fatol=1e-12))
    hw_p, hg = float(res.x[0]), float(res.x[1])
    jump = abs(_drude(lam_hi, eps_inf, hw_p, hg) - fn(lam_hi)) / abs(fn(lam_hi))
    return hw_p, hg, float(jump)


class AuModel:
    """Золото на всём диапазоне: таблица в своих пределах, Друде за их верхней границей.

    ir = 'fit': Друде с eps_inf = 1, hw_p и hgamma подобраны по последним 500 нм таблицы;
    ir = 'ordal': Друде по Ordal 1985 (hw_p = 9.03 эВ, hgamma = 0.0267 эВ) без подгонки - для оценки
    чувствительности результата к ИК-модели золота.
    """

    def __init__(self, table: str = "lerer", ir: str = "fit"):
        self.table = table
        self.ir = ir
        if table == "lerer":
            self.fn, (self.lo, self.hi) = eps_au_lerer, LERER_RANGE
        elif table == "jc":
            self.fn, (self.lo, self.hi) = eps_au_jc, JC_RANGE
        else:
            raise ValueError(table)
        self.eps_inf = 1.0
        if ir == "fit":
            self.hw_p, self.gamma, self.jump = _drude_fit(self.fn, self.hi - 500.0, self.hi)
        elif ir == "ordal":
            self.hw_p, self.gamma = HW_P_AU_EV, HGAMMA_ORDAL_EV
            self.jump = abs(_drude(self.hi, 1.0, self.hw_p, self.gamma) - self.fn(self.hi)) / abs(self.fn(self.hi))
        else:
            raise ValueError(ir)

    def __call__(self, lam_nm: float) -> complex:
        if lam_nm <= self.hi:
            if lam_nm < self.lo:
                raise ValueError("lambda ниже таблицы: %.1f нм" % lam_nm)
            return self.fn(lam_nm)
        return _drude(lam_nm, self.eps_inf, self.hw_p, self.gamma)

    def describe(self) -> str:
        return ("Au: таблица %s до %.0f нм, дальше Друде (%s) hw_p=%.3f эВ, hgamma=%.4f эВ, eps_inf=%.1f; "
                "скачок на стыке %.1f %%" % (self.table, self.hi, self.ir, self.hw_p, self.gamma, self.eps_inf,
                                              100 * self.jump))


# ---------------------------------------------------------------------------
# Правила смешивания (оптическая конвенция)
# ---------------------------------------------------------------------------

def k_contrast(eps_h: complex, eps_p: complex) -> complex:
    return (eps_p - eps_h) / (eps_p + 2.0 * eps_h)


def mg(eps_h: complex, eps_p: complex, f: float) -> complex:
    """Полный Максвелл-Гарнетт; в COMPOSIT.c назван формулой Клаузиуса-Моссотти."""
    q = f * k_contrast(eps_h, eps_p)
    return eps_h * (1.0 + 2.0 * q) / (1.0 - q)


def mg_multi(eps_h: complex, populations: list[tuple[complex, float]]) -> complex:
    """Максвелл-Гарнетт для нескольких сортов частиц - исправленная COMPOSITE_3."""
    q = sum(f * k_contrast(eps_h, eps_p) for eps_p, f in populations)
    return eps_h * (1.0 + 2.0 * q) / (1.0 - q)


def dilute(eps_h: complex, eps_p: complex, f: float) -> complex:
    return eps_h + 3.0 * f * eps_h * k_contrast(eps_h, eps_p)


def lerer_formula(eps_h: complex, eps_p: complex, f: float) -> complex:
    """Запись Лерера 2026-08-19: 1 + (1-C)(eps_h-1) + C chi2, chi2 = 3 eps_h K."""
    return 1.0 + (1.0 - f) * (eps_h - 1.0) + 3.0 * f * eps_h * k_contrast(eps_h, eps_p)


def bruggeman(eps_h: complex, eps_p: complex, f: float) -> complex:
    """Симметричный Бруггеман; из двух корней берётся пассивный (Im >= 0), ближайший к eps_h."""
    b = (3.0 * f - 1.0) * eps_p + (2.0 - 3.0 * f) * eps_h
    disc = cmath.sqrt(b * b + 8.0 * eps_p * eps_h)
    roots = [(b + disc) / 4.0, (b - disc) / 4.0]
    passive = [r for r in roots if r.imag >= -1e-12]
    cand = passive if passive else roots
    return min(cand, key=lambda z: abs(z - eps_h))


def lorentz_lorenz(eps_h: complex, eps_p: complex, f: float) -> complex:
    s = (1.0 - f) * (eps_h - 1.0) / (eps_h + 2.0) + f * (eps_p - 1.0) / (eps_p + 2.0)
    return (1.0 + 2.0 * s) / (1.0 - s)


def mlwa_mg(eps_h: complex, eps_p: complex, f: float, lam_nm: float, radius_nm: float) -> complex:
    """Максвелл-Гарнетт через поляризуемость сферы с поправкой MLWA (динамическая деполяризация
    и реакция излучения), как в compare_effective_models.py; здесь сразу в оптической конвенции."""
    r = radius_nm * 1e-9
    n_h = cmath.sqrt(eps_h).real
    k_h = 2.0 * math.pi * n_h / (lam_nm * 1e-9)
    beta0 = r ** 3 * k_contrast(eps_h, eps_p)
    beta = beta0 / (1.0 - (k_h ** 2 / r) * beta0 - 1j * (2.0 / 3.0) * k_h ** 3 * beta0)
    eta = f * beta / r ** 3
    return eps_h * (1.0 + 2.0 * eta) / (1.0 - eta)


MIXING = {
    "MG": lambda eh, ep, f, lam: mg(eh, ep, f),
    "Bruggeman": lambda eh, ep, f, lam: bruggeman(eh, ep, f),
    "MLWA15": lambda eh, ep, f, lam: mlwa_mg(eh, ep, f, lam, 15.0),
    "Lerer": lambda eh, ep, f, lam: lerer_formula(eh, ep, f),
    "dilute": lambda eh, ep, f, lam: dilute(eh, ep, f),
    "LL": lambda eh, ep, f, lam: lorentz_lorenz(eh, ep, f),
}
MIXING_RU = {
    "MG": "Максвелл-Гарнетт (у Лерера - Клаузиус-Моссотти)",
    "Bruggeman": "Бруггеман",
    "MLWA15": "Максвелл-Гарнетт с поправкой MLWA, R = 15 нм",
    "Lerer": "формула Лерера",
    "dilute": "первый порядок по концентрации",
    "LL": "Лоренц-Лоренц",
}


def n_of(eps: complex) -> complex:
    n = cmath.sqrt(eps)
    return n if n.imag >= 0 else -n


# ---------------------------------------------------------------------------
# Описание стопки слоёв и сборка объекта grcwa
# ---------------------------------------------------------------------------

@dataclass
class Uniform:
    t_um: float
    eps: complex
    name: str = ""


@dataclass
class Pillar:
    """Слой с решёткой столбиков одного сечения в среде eps_out.

    shape: 'square' | 'circle' | 'ellipse'; ax, ay - полные размеры сечения по x и y, мкм
    (для квадрата и круга ax = ay). angle_deg - поворот эллипса относительно оси x.
    """
    t_um: float
    ax_um: float
    ay_um: float
    eps_in: complex
    eps_out: complex = 1.0 + 0j
    shape: str = "square"
    angle_deg: float = 0.0
    name: str = ""


def cone_layers(t_um: float, a_bottom_um: float, a_top_um: float, eps_in: complex, n_slices: int,
                shape: str = "square", eps_out: complex = 1.0 + 0j, name: str = "cone") -> list[Pillar]:
    """Усечённый конус как лестница из n_slices слоёв; размер каждой ступеньки берётся в её середине.
    Список идёт сверху вниз (первым - самый узкий слой), как и вся стопка в grcwa."""
    out = []
    for i in range(n_slices):
        zc = (i + 0.5) / n_slices          # 0 - верх, 1 - низ
        a = a_top_um + (a_bottom_um - a_top_um) * zc
        out.append(Pillar(t_um / n_slices, a, a, eps_in, eps_out, shape, 0.0, "%s%d" % (name, i)))
    return out


def pillar_mask(p: Pillar, period_um: float, nx: int) -> np.ndarray:
    x = (np.arange(nx) + 0.5) / nx * period_um - period_um / 2
    xx, yy = np.meshgrid(x, x, indexing="ij")
    if p.shape == "square":
        return (np.abs(xx) <= p.ax_um / 2) & (np.abs(yy) <= p.ay_um / 2)
    if p.shape in ("circle", "ellipse"):
        th = math.radians(p.angle_deg)
        xr = xx * math.cos(th) + yy * math.sin(th)
        yr = -xx * math.sin(th) + yy * math.cos(th)
        return (xr / (p.ax_um / 2)) ** 2 + (yr / (p.ay_um / 2)) ** 2 <= 1.0
    raise ValueError(p.shape)


def fill_fraction(p: Pillar, period_um: float, nx: int = 400) -> float:
    return float(pillar_mask(p, period_um, nx).mean())


@dataclass
class Solution:
    R: float
    T: float
    A: float                       # 1 - R - T: поглощение в конечных слоях
    A_layers: list[float]          # поглощение по слоям стопки (без полупространств)
    obj: object = field(repr=False, default=None)


def build(lam_nm: float, period_um: float, layers: list, eps_super: complex = 1.0 + 0j,
          eps_sub: complex = N_SUB ** 2 + 0j, nG: int = 121, pol: str = "x", nx: int = 120,
          theta_deg: float = 0.0, phi_deg: float = 0.0):
    """Собирает объект grcwa для стопки layers (сверху вниз) между полупространствами.

    pol = 'x' - E вдоль x (p_amp = 1), 'y' - E вдоль y (s_amp = 1); при нормальном падении на
    квадратную решётку оба совпадают. theta_deg - угол падения; строго нулевой угол делает
    матрицу kp вырожденной, поэтому подставляется 1e-4 рад.
    """
    lam = lam_nm / 1000.0
    theta = math.radians(theta_deg) if theta_deg > 0 else 1e-4
    obj = grcwa.obj(nG, [period_um, 0.0], [0.0, period_um], 1.0 / lam, theta, math.radians(phi_deg), verbose=0)
    obj.Add_LayerUniform(1.0, eps_super)
    grids = []
    for L in layers:
        if isinstance(L, Uniform):
            obj.Add_LayerUniform(L.t_um, L.eps)
        else:
            obj.Add_LayerGrid(L.t_um, nx, nx)
    obj.Add_LayerUniform(1.0, eps_sub)
    obj.Init_Setup()
    for L in layers:
        if isinstance(L, Pillar):
            m = pillar_mask(L, period_um, nx)
            grids.append(np.where(m, L.eps_in, L.eps_out).astype(complex))
    if grids:
        obj.GridLayer_geteps(np.concatenate([g.flatten() for g in grids]))
    if pol == "x":
        obj.MakeExcitationPlanewave(1.0, 0.0, 0.0, 0.0, order=0)
    else:
        obj.MakeExcitationPlanewave(0.0, 0.0, 1.0, 0.0, order=0)
    return obj, grids


ABS_FACTOR = 1.0   # калибруется в self_test: доля поглощённой мощности = ABS_FACTOR * omega * Volume_integral


def layer_absorption_volume(obj, layer_index: int, im_eps, nx: int = 120) -> float:
    """Поглощение в слое через объёмный интеграл Im eps |E|^2 решателя (Volume_integral).

    Для однородных слоёв совпадает с разностью потоков точно; для слоя с решёткой правило Лорана
    для произведения разрывных Im eps(x, y) и E_x, E_y занижает интеграл (в самопроверке - на 27...32 %
    при 121 гармонике), поэтому в расчётах используется layer_absorption() по потокам.
    """
    nG = obj.nG
    if np.isscalar(im_eps):
        M = float(im_eps) * np.eye(nG, dtype=complex)
    else:
        M = get_conv(1.0 / (nx * nx), np.asarray(im_eps, dtype=complex), obj.G)
    val = obj.Volume_integral(layer_index, M, M, M, normalize=1)
    return ABS_FACTOR * float(np.real(obj.omega * val))


def layer_absorption(obj, layer_index: int, im_eps=None, nx: int = 120) -> float:
    """Доля падающей мощности, поглощённая в слое layer_index (0 = верхнее полупространство):
    разность нормальных потоков Пойнтинга через верхнюю и нижнюю границы слоя, посчитанных из тех же
    амплитуд, что дают R и T. im_eps не используется и оставлен для совместимости вызова."""
    q, kp, phi = obj.q_list[layer_index], obj.kp_list[layer_index], obj.phi_list[layer_index]
    t = obj.thickness_list[layer_index]
    a_top, b_top = obj.GetAmplitudes(layer_index, 0.0)
    f1, b1 = _rc.GetZPoyntingFlux(a_top, b_top, obj.omega, kp, phi, q)
    a_bot, b_bot = obj.GetAmplitudes(layer_index, t)
    f2, b2 = _rc.GetZPoyntingFlux(a_bot, b_bot, obj.omega, kp, phi, q)
    net_top = np.real(f1 + b1)
    net_bot = np.real(f2 + b2)
    return float(np.real((net_top - net_bot) * obj.normalization))


def solve(lam_nm: float, period_um: float, layers: list, want_layers: bool = False, **kw) -> Solution:
    obj, grids = build(lam_nm, period_um, layers, **kw)
    R, T = obj.RT_Solve(normalize=1)
    R, T = float(np.real(R)), float(np.real(T))
    A_layers = []
    if want_layers:
        gi = 0
        nx = kw.get("nx", 120)
        for i, L in enumerate(layers):
            if isinstance(L, Uniform):
                A_layers.append(layer_absorption(obj, i + 1, L.eps.imag, nx))
            else:
                A_layers.append(layer_absorption(obj, i + 1, grids[gi].imag, nx))
                gi += 1
    return Solution(R, T, 1.0 - R - T, A_layers, obj)


def mean_e2_ratio(A_layer: float, lam_nm: float, im_eps: float, t_um: float, fill: float) -> float:
    """<|E|^2> в области композита слоя, отнесённое к |E0|^2 падающей волны в воздухе.

    A = k0 Im eps <|E|^2> t fill / |E0|^2  (с = eps0 = 1), откуда и берётся отношение."""
    k0 = 2.0 * math.pi / (lam_nm / 1000.0)
    return A_layer / (k0 * im_eps * t_um * fill)


# ---------------------------------------------------------------------------
# Стандартные конструкции
# ---------------------------------------------------------------------------

def stack_pillar_absorber(lam_nm: float, au, period_nm: float = 400.0, pillar_nm: float = 250.0,
                          h_pil_nm: float = 100.0, layers_comp=((0.10, 2000.0),), spacer_nm: float = 325.0,
                          mirror: bool = False, shape: str = "square", pillar_layers=None,
                          eps_super: complex = 1.0 + 0j, mixing=mg):
    """Столбиковый поглотитель Лерера (описание конструкции 21.08.2026) и его переделки.

    layers_comp - слои композита сверху вниз (доля Au, толщина нм); spacer - прослойка n = 1.45;
    mirror - золото вместо подложки. pillar_layers - готовый список слоёв столбиков (например конус)
    вместо одного слоя; их композит берётся по первой доле из layers_comp.
    Возвращает (layers, eps_sub)."""
    eps_h = complex(N_HOST ** 2, 0.0)
    eps_c_top = mixing(eps_h, au(lam_nm), layers_comp[0][0])
    layers = []
    if pillar_layers is not None:
        for p in pillar_layers:
            layers.append(Pillar(p.t_um, p.ax_um, p.ay_um, eps_c_top, eps_super, p.shape, p.angle_deg, p.name))
    elif h_pil_nm > 0:
        layers.append(Pillar(h_pil_nm / 1000.0, pillar_nm / 1000.0, pillar_nm / 1000.0, eps_c_top, eps_super,
                             shape, 0.0, "pillars"))
    for f, t in layers_comp:
        layers.append(Uniform(t / 1000.0, mixing(eps_h, au(lam_nm), f), "comp%.0f" % (100 * f)))
    if spacer_nm > 0:
        layers.append(Uniform(spacer_nm / 1000.0, N_SUB ** 2 + 0j, "spacer"))
    eps_sub = au(lam_nm) if mirror else N_SUB ** 2 + 0j
    return layers, eps_sub


ORIGINAL = dict(period_nm=400.0, pillar_nm=250.0, h_pil_nm=100.0, layers_comp=((0.10, 2000.0),),
                spacer_nm=325.0, mirror=False)
IMPROVED = dict(period_nm=200.0, pillar_nm=125.0, h_pil_nm=100.0, layers_comp=((0.10, 500.0), (0.20, 1500.0)),
                spacer_nm=50.0, mirror=True)
IMPROVED_NO_MIRROR = dict(period_nm=200.0, pillar_nm=125.0, h_pil_nm=100.0,
                          layers_comp=((0.10, 500.0), (0.20, 1500.0)), spacer_nm=325.0, mirror=False)


def stack_cylinder_metasurface(lam_nm: float, eps_comp_cyl: complex, eps_comp_film: complex,
                               period_nm: float = 500.0, diam_nm: float = 400.0, h_cyl_nm: float = 100.0,
                               h_film_nm: float = 100.0, spacer_nm: float = 0.0):
    """Метаповерхность статьи Лерер АМ_1: цилиндры композита в воздухе над слоем композита на подложке 1.45.
    spacer_nm > 0 добавляет прослойку n = 1.45 (входной файл письма 19.08.2026: 325 нм)."""
    layers = [Pillar(h_cyl_nm / 1000.0, diam_nm / 1000.0, diam_nm / 1000.0, eps_comp_cyl, 1.0 + 0j, "circle",
                     0.0, "cyl"),
              Uniform(h_film_nm / 1000.0, eps_comp_film, "film")]
    if spacer_nm > 0:
        layers.append(Uniform(spacer_nm / 1000.0, N_SUB ** 2 + 0j, "spacer"))
    return layers, N_SUB ** 2 + 0j


# ---------------------------------------------------------------------------
# Аналитическая матрица переноса для проверки однородных пределов
# ---------------------------------------------------------------------------

def tmm_rta(lam_nm: float, layers_n_t: list[tuple[complex, float]], n_in: complex, n_out: complex):
    """R, T, A стопки однородных слоёв [(n, толщина нм)] при нормальном падении, из среды n_in в n_out."""
    m = np.eye(2, dtype=complex)
    for n, t in layers_n_t:
        d = 2.0 * math.pi * n * t / lam_nm
        c, s = cmath.cos(d), cmath.sin(d)
        m = m @ np.array([[c, -1j * s / n], [-1j * n * s, c]], dtype=complex)
    a = (m[0, 0] + m[0, 1] * n_out) * n_in
    b = m[1, 0] + m[1, 1] * n_out
    r = (a - b) / (a + b)
    t = 2.0 * n_in / (a + b)
    R = abs(r) ** 2
    T = (complex(n_out).real / complex(n_in).real) * abs(t) ** 2
    return R, T, 1.0 - R - T


def self_test(verbose: bool = True) -> float:
    """Калибровка множителя поглощения по слоям на однородной стопке и проверка суммы по слоям."""
    global ABS_FACTOR
    au = AuModel("lerer")
    lam = 550.0
    eps_h = complex(N_HOST ** 2, 0.0)
    ec = mg(eps_h, au(lam), 0.10)
    layers = [Uniform(0.10, ec, "c1"), Uniform(0.325, N_SUB ** 2 + 0j, "sp"), Uniform(0.05, ec, "c2")]
    ABS_FACTOR = 1.0
    sol = solve(lam, 0.4, layers, want_layers=True, nG=25)
    R_t, T_t, A_t = tmm_rta(lam, [(n_of(ec), 100.0), (N_SUB, 325.0), (n_of(ec), 50.0)], 1.0, N_SUB)
    vol = [layer_absorption_volume(sol.obj, i + 1, L.eps.imag) for i, L in enumerate(layers)]
    factor = (1.0 - sol.R - sol.T) / sum(vol)
    ABS_FACTOR = round(factor, 6)
    if verbose:
        print("[self_test] однородная стопка при 550 нм: grcwa R=%.6f T=%.6f A=%.6f | матрица переноса "
              "R=%.6f T=%.6f A=%.6f" % (sol.R, sol.T, sol.A, R_t, T_t, A_t))
        print("[self_test] сумма omega*Volume_integral по слоям = %.6f при A = %.6f -> множитель %.6f"
              % (sum(vol), sol.A, factor))
        print("[self_test] по потокам: сумма по слоям = %.6f, по слоям %s"
              % (sum(sol.A_layers), " ".join("%.5f" % a for a in sol.A_layers)))
    assert abs(sol.R - R_t) < 2e-4 and abs(sol.T - T_t) < 2e-4, "grcwa не совпал с матрицей переноса"
    assert abs(sum(sol.A_layers) - sol.A) < 1e-6, "разность потоков не сходится с 1 - R - T"
    return factor


if __name__ == "__main__":
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    f = self_test()
    print("ABS_FACTOR =", ABS_FACTOR)
    au = AuModel("lerer")
    print(au.describe())
    print(AuModel("jc").describe())
    print(AuModel("lerer", "ordal").describe())
    for lam in (1900.0, 2000.0, 2100.0, 3000.0, 10000.0):
        print("  Au Lerer/Drude %6.0f нм: %s" % (lam, au(lam)))
    # столбики: доля заполнения и проверка суммы поглощения по слоям на решётке
    eps_h = complex(N_HOST ** 2, 0.0)
    lam = 550.0
    layers, eps_sub = stack_pillar_absorber(lam, au, **ORIGINAL)
    sol = solve(lam, 0.4, layers, eps_sub=eps_sub, want_layers=True, nG=121)
    print("[решётка] A = %.5f, сумма по слоям = %.5f, по слоям: %s"
          % (sol.A, sum(sol.A_layers), " ".join("%.4f" % a for a in sol.A_layers)))
    print("[решётка] заполнение квадрата 250/400: %.4f" % fill_fraction(layers[0], 0.4))
