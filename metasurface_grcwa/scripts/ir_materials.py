# -*- coding: utf-8 -*-
"""Оптические постоянные наночастиц для статьи об ИК-диапазоне (400-2000 нм), конвенция exp(-i w t), Im eps >= 0.

Основные наборы - из публикаций, таблицы программы Лерера - для сравнения (решение 05.11.2026):
  TiN_R17   Reddy et al., arXiv:1702.03053, табл. S1: плёнка 200 нм на сапфире, 23 C, Друде + 2 Лоренца;
  TiN_S22   Saha et al., arXiv:2203.14170, табл. S1: плёнка 200 нм на Si, Друде + 2 Лоренца;
  TiN_P84   Pflüger et al. 1984 из composite_ema: в ИК всего узлы 826, 1240, 2480 нм - только для оценки разброса;
  ZrN_N11   Naik et al. 2011, рис. 4, оцифровка из composite_ema (плёнка на сапфире): 300-2000 нм;
  AZO_S22   Saha et al., табл. 2 приложения: плёнка 63 нм, мишень 2 % Al, Друде + 1 Лоренц;
  AZO_S22b  то же, плёнка 57 нм: оценка разброса;
  Al_R98    Rakić et al. 1998, модель Лоренца-Друде, табл. 2;
  Cu_JC     Johnson, Christy 1972 до 1937 нм, далее Друде, подобранный по последним 500 нм таблицы;
  Cu_R98    Rakić et al. 1998, модель Лоренца-Друде;
  Al_L, Cu_L  таблицы программы Лерера (Al_.c, libmat 1).

Формулы моделей: Reddy - уравнение (1) статьи; Saha - уравнение (S1), слагаемое A_n / (E_n^2 - E^2 - i B_n E).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "vibr_t" / "scripts"))
import lerer_grcwa as L  # noqa: E402
import vibr_materials as VM  # noqa: E402

ce = L.ce
HC = 1239.841984


def _e(lam_nm: float) -> float:
    return HC / lam_nm


def drude_lorentz(lam_nm: float, eps_inf: float, wp2: float, gamma: float, osc) -> complex:
    """eps = eps_inf - wp2 / (E^2 + i gamma E) + sum A / (E0^2 - E^2 - i G E), энергии в эВ."""
    E = _e(lam_nm)
    z = eps_inf - wp2 / (E * E + 1j * gamma * E)
    for A, G, E0 in osc:
        z += A / (E0 * E0 - E * E - 1j * G * E)
    return z


def rakic_ld(lam_nm: float, wp: float, f0: float, g0: float, osc) -> complex:
    """Rakić et al. 1998: eps = 1 - f0 wp^2 / (E (E + i g0)) + sum f wp^2 / (Ej^2 - E^2 - i E Gj); osc = [(f, G, Ej)]."""
    E = _e(lam_nm)
    z = 1.0 - f0 * wp * wp / (E * (E + 1j * g0))
    for f, G, Ej in osc:
        z += f * wp * wp / (Ej * Ej - E * E - 1j * E * G)
    return z


def eps_tin_reddy(lam):
    return drude_lorentz(lam, 5.18, 7.38 ** 2, 0.26, [(42.2, 1.42, 4.07), (2.26, 0.87, 2.02)])


def eps_tin_saha(lam):
    return drude_lorentz(lam, 3.41, 40.152, 0.227, [(28.517, 1.003, 4.117), (391.66, 99.451, 4.748)])


def eps_azo_saha63(lam):
    return drude_lorentz(lam, 3.224, 2.64, 0.135, [(7.412, 0.019, 4.258)])


def eps_azo_saha57(lam):
    return drude_lorentz(lam, 2.392, 2.519, 0.146, [(21.544, 0.006, 4.386)])


AL_LD = (14.98, 0.523, 0.047, [(0.227, 0.333, 0.162), (0.050, 0.312, 1.544), (0.166, 1.351, 1.808), (0.030, 3.382, 3.473)])
CU_LD = (10.83, 0.575, 0.030, [(0.061, 0.378, 0.291), (0.104, 1.056, 2.957), (0.723, 3.213, 5.300), (0.638, 4.305, 11.18)])


def eps_al_rakic(lam):
    return rakic_ld(lam, *AL_LD)


def eps_cu_rakic(lam):
    return rakic_ld(lam, *CU_LD)


_CU_HI = ce.MATERIAL_RANGE_NM[1][1]
_CU_DRUDE = L._drude_fit(L.eps_cu_jc, _CU_HI - 500.0, _CU_HI)


def eps_cu_jc(lam):
    if lam <= _CU_HI:
        return L.eps_cu_jc(lam)
    return L._drude(lam, 1.0, _CU_DRUDE[0], _CU_DRUDE[1])


def _load(name):
    a = np.loadtxt(L.DATA / name, comments="#")
    return a[:, 0], a[:, 1], a[:, 2]


_AL_L = _load("Al_lerer_nk.txt")


def eps_al_lerer(lam):
    return complex(float(np.interp(lam, _AL_L[0], _AL_L[1])), float(np.interp(lam, _AL_L[0], _AL_L[2]))) ** 2


def eps_cu_lerer(lam):
    return VM.libmat(1, float(lam)).conjugate()


def eps_tin_pfluger(lam):
    return ce.eps_tin(float(lam)).conjugate()


def eps_zrn_naik(lam):
    return ce.eps_zrn(float(lam)).conjugate()


MATERIALS = {
    "TiN_R17": eps_tin_reddy, "TiN_S22": eps_tin_saha, "TiN_P84": eps_tin_pfluger, "ZrN_N11": eps_zrn_naik,
    "AZO_S22": eps_azo_saha63, "AZO_S22b": eps_azo_saha57, "Al_R98": eps_al_rakic, "Al_L": eps_al_lerer,
    "Cu_JC": eps_cu_jc, "Cu_R98": eps_cu_rakic, "Cu_L": eps_cu_lerer,
}
MAIN = ("TiN_R17", "ZrN_N11", "AZO_S22", "Al_R98", "Cu_JC")
LABEL = {"TiN_R17": "TiN (Reddy 2017)", "TiN_S22": "TiN (Saha 2022)", "TiN_P84": "TiN (Pflüger 1984)", "ZrN_N11": "ZrN (Naik 2011)",
         "AZO_S22": "AZO (Saha 2022, 63 нм)", "AZO_S22b": "AZO (Saha 2022, 57 нм)", "Al_R98": "Al (Rakić 1998)", "Al_L": "Al (программа Лерера)",
         "Cu_JC": "Cu (Johnson, Christy)", "Cu_R98": "Cu (Rakić 1998)", "Cu_L": "Cu (программа Лерера)"}


def frohlich(name, eps_h=L.N_HOST ** 2, lo=300.0, hi=4000.0):
    """Длина волны, где Re eps_p = -2 eps_h (нм), или None; вне таблиц ZrN (2000 нм) и Pflüger (2480 нм) не считается."""
    fn = MATERIALS[name]
    if name == "ZrN_N11":
        hi = 2000.0
    if name == "TiN_P84":
        hi = 2479.0
    if name in ("Cu_L", "Al_L"):
        hi = 2000.0
    lams = np.arange(lo, hi, 1.0)
    d = np.array([fn(float(x)).real + 2.0 * eps_h for x in lams])
    i = np.where(np.sign(d[:-1]) * np.sign(d[1:]) < 0)[0]
    return [float(lams[j]) for j in i]


def composite(names, fracs, lam, eps_h=complex(L.N_HOST ** 2, 0.0)):
    """Максвелл-Гарнетт для одного или нескольких сортов частиц."""
    return L.mg_multi(eps_h, [(MATERIALS[n](lam), f) for n, f in zip(names, fracs)])


def _check():
    """Контрольные значения модели Reddy и Rakić по формуле, вручную при 1550 нм."""
    E = HC / 1550.0
    drude = -7.38 ** 2 / (E * E + 1j * 0.26 * E)
    assert abs(eps_tin_reddy(1550.0) - (5.18 + drude + 42.2 / (4.07 ** 2 - E * E - 1j * 1.42 * E) + 2.26 / (2.02 ** 2 - E * E - 1j * 0.87 * E))) < 1e-9
    for n, fn in MATERIALS.items():
        z = fn(1000.0)
        assert z.imag >= 0, n


if __name__ == "__main__":
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    _check()
    print("Cu: Друде за %g нм: hw_p=%.3f эВ, hg=%.4f эВ, скачок на стыке %.1f %%" % (_CU_HI, _CU_DRUDE[0], _CU_DRUDE[1], 100 * _CU_DRUDE[2]))
    eh = L.N_HOST ** 2
    print("Основа n = %.2f, eps_h = %.3f; условие Фрёлиха: Re eps_p = %.3f; f = 0.1" % (L.N_HOST, eh, -2 * eh))
    lams = (500.0, 800.0, 1000.0, 1550.0, 2000.0)
    print("%-10s" % "материал" + "".join("%18s" % ("eps @%g" % x) for x in lams) + "   Фрёлих, нм")
    for n, fn in MATERIALS.items():
        print("%-10s" % n + "".join("%18s" % ("%.2f%+.2fi" % (fn(x).real, fn(x).imag)) for x in lams) + "   " + str(frohlich(n)))
    print()
    print("Im eps_eff и ослабление за проход слоя 2 мкм, дБ (f = 0.1, Максвелл-Гарнетт)")
    print("%-10s" % "материал" + "".join("%18s" % ("%g нм" % x) for x in lams))
    for n in MATERIALS:
        row = ""
        for x in lams:
            ee = composite([n], [0.1], x)
            nn = L.n_of(ee)
            att = 4.343 * 4 * 3.14159265 * nn.imag / x * 2000.0
            row += "%18s" % ("%.3f / %.1f дБ" % (ee.imag, att))
        print("%-10s" % n + row)
