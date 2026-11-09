# -*- coding: utf-8 -*-
"""Общее для расчётов статьи об ИК-диапазоне: составы композита, спектры, показатели полосы, папка данных."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import ir_materials as M  # noqa: E402
import lerer_grcwa as L  # noqa: E402

# состав: {имя: ((материалы), (объёмные доли))}; основа n = 1.77, формула КМ (Максвелл-Гарнетт)
COMBOS = {
    "TiN": (("TiN_R17",), (0.10,)),
    "ZrN": (("ZrN_N11",), (0.10,)),
    "AZO": (("AZO_S22",), (0.10,)),
    "Al": (("Al_R98",), (0.10,)),
    "Cu": (("Cu_JC",), (0.10,)),
    "Cu+AZO_05": (("Cu_JC", "AZO_S22"), (0.05, 0.05)),
    "TiN+AZO_05": (("TiN_R17", "AZO_S22"), (0.05, 0.05)),
    "ZrN+AZO_05": (("ZrN_N11", "AZO_S22"), (0.05, 0.05)),
    "Cu+TiN_05": (("Cu_JC", "TiN_R17"), (0.05, 0.05)),
    "Cu+AZO_10": (("Cu_JC", "AZO_S22"), (0.10, 0.10)),
    "TiN+AZO_10": (("TiN_R17", "AZO_S22"), (0.10, 0.10)),
    "ZrN+AZO_10": (("ZrN_N11", "AZO_S22"), (0.10, 0.10)),
    # разброс данных и таблицы программы Лерера: те же составы с другими наборами
    "TiN_S22": (("TiN_S22",), (0.10,)),
    "TiN_P84": (("TiN_P84",), (0.10,)),
    "AZO_S22b": (("AZO_S22b",), (0.10,)),
    "Al_L": (("Al_L",), (0.10,)),
    "Cu_L": (("Cu_L",), (0.10,)),
    "Cu_R98": (("Cu_R98",), (0.10,)),
    "TiN_S22+AZO_10": (("TiN_S22", "AZO_S22"), (0.10, 0.10)),
    "TiN_P84+AZO_10": (("TiN_P84", "AZO_S22"), (0.10, 0.10)),
    "TiN+AZO57_10": (("TiN_R17", "AZO_S22b"), (0.10, 0.10)),
    "ZrN+AZO57_10": (("ZrN_N11", "AZO_S22b"), (0.10, 0.10)),
}
LAMS = np.arange(400.0, 2001.0, 50.0)
EPS_SUB = L.N_SUB ** 2 + 0j


def eps_combo(name, lam):
    names, fracs = COMBOS[name]
    return M.composite(names, fracs, float(lam))


def data_dir():
    """Папка данных статьи в Obsidian-служебное (имя корня Диска отличается на разных ПК)."""
    root = next(p for p in Path.home().glob("Мой диск*/Obsidian-служебное"))
    d = root / "PhD/Scientific_study/Scientific_study/paper_analysis/from_Lerer/от Лерера/Статья 2 ИК"
    d.mkdir(parents=True, exist_ok=True)
    return d


def stack(name, lam, h_layer_nm, ax_nm=0.0, h_pillar_nm=0.0, buffer_nm=325.0):
    """Столбики (если есть) + слой композита + буфер n = 1.45; подложка n = 1.45."""
    e = eps_combo(name, lam)
    layers = []
    if h_pillar_nm > 0 and ax_nm > 0:
        layers.append(L.Pillar(h_pillar_nm / 1000.0, ax_nm / 1000.0, ax_nm / 1000.0, e, 1.0 + 0j, "square", 0.0, "pillars"))
    layers.append(L.Uniform(h_layer_nm / 1000.0, e, "composite"))
    if buffer_nm > 0:
        layers.append(L.Uniform(buffer_nm / 1000.0, EPS_SUB, "buffer"))
    return layers


def solve_point(name, lam, period_nm, h_layer_nm, ax_nm, h_pillar_nm, nG=81, nx=120):
    return L.solve(float(lam), period_nm / 1000.0, stack(name, lam, h_layer_nm, ax_nm, h_pillar_nm), eps_sub=EPS_SUB, nG=nG, nx=nx)


def tmm_spectrum(name, h_layer_nm, lams=LAMS):
    out = []
    for lam in lams:
        n = L.n_of(eps_combo(name, lam))
        out.append(L.tmm_rta(float(lam), [(n, h_layer_nm), (L.N_SUB, 325.0)], 1.0, L.N_SUB))
    return np.array(out)  # R, T, A


def band(lams, A, level=0.9):
    """Наибольший непрерывный участок A >= level: (ширина нм, начало, конец); границы - линейной интерполяцией."""
    lams = np.asarray(lams, float)
    A = np.asarray(A, float)
    best = (0.0, None, None)
    i, n = 0, len(lams)
    while i < n:
        if A[i] < level:
            i += 1
            continue
        j = i
        while j + 1 < n and A[j + 1] >= level:
            j += 1
        lo = lams[i] if i == 0 else lams[i - 1] + (level - A[i - 1]) / (A[i] - A[i - 1]) * (lams[i] - lams[i - 1])
        hi = lams[j] if j == n - 1 else lams[j] + (level - A[j]) / (A[j + 1] - A[j]) * (lams[j + 1] - lams[j])
        if hi - lo > best[0]:
            best = (float(hi - lo), float(lo), float(hi))
        i = j + 1
    return best
